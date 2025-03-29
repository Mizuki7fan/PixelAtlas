import os
import subprocess
import shutil
import json

def fix_compile_commands(build_dir: str):
    """
    修复 compile_commands.json 中的 MSVC 参数，使其兼容 clangd
    """
    compile_commands_path = os.path.join(build_dir, "compile_commands.json")
    
    if not os.path.exists(compile_commands_path):
        print(f"Error: {compile_commands_path} not found")
        return

    # 读取并修复 JSON 内容
    with open(compile_commands_path, "r", encoding="utf-8") as f:
        entries = json.load(f)

    for entry in entries:
        cmd = entry["command"]
        # 替换 /external:I 为 -isystem
        cmd = cmd.replace("/external:I", "-isystem")
        # 替换路径分隔符（Windows 反斜杠 → 正斜杠）
        cmd = cmd.replace("\\", "/")
        entry["command"] = cmd

    # 写回文件
    with open(compile_commands_path, "w", encoding="utf-8") as f:
        json.dump(entries, f, indent=2, ensure_ascii=False)
    
    print(f"Fixed: {compile_commands_path}")

def get_project_root(current_file: str) -> str:
    """根据脚本路径动态推导项目根目录（脚本位于 /pa-script/etc/）"""
    script_dir = os.path.dirname(os.path.abspath(current_file))  # G:/PixelAtlas/pa-script/etc
    project_root = os.path.abspath(os.path.join(script_dir, "../../"))  # 上溯两级到项目根目录
    return project_root

def generate_compile_commands(
    source_dir: str, 
    build_dir: str,
    generator: str,
    vs_tools_path: str
):
    """
    快速生成 compile_commands.json（无需构建项目）
    
    Args:
        source_dir: 项目源码目录（包含 CMakeLists.txt）
        build_dir: 构建目录（默认在源码目录下创建 build）
        generator: CMake 生成器（推荐 Ninja）
        vs_tools_path: Visual Studio 环境脚本路径（仅 Windows 需要）
    """
    # 清理并创建构建目录
    build_path = os.path.join(source_dir, build_dir)
    if os.path.exists(build_path):
        shutil.rmtree(build_path)
    os.makedirs(build_path, exist_ok=True)

    # 配置 CMake 命令
    cmake_cmd = [
        "cmake",
        f"-G{generator}",
        "--log-level=WARNING", # 仅显示警告和错误
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-B", build_path,
        "-S", source_dir
    ]

    # Windows 下需要加载 MSVC 环境
    if os.name == "nt" and generator == "Ninja":
        # 调用 vcvars 设置环境变量
        subprocess.run(f'"{vs_tools_path}" && cmake {" ".join(cmake_cmd[1:])}', shell=True, cwd=build_path)
    else:
        subprocess.run(cmake_cmd, cwd=build_path)

    # 检查生成结果
    compile_commands = os.path.join(build_path, "compile_commands.json")
    if os.path.exists(compile_commands):
        print(f"Success: {compile_commands}")
        #fix_compile_commands(build_path)
    else:
        print("Error: Failed to generate compile_commands.json")

if __name__ == "__main__":
    project_root = get_project_root(__file__)
    build_dir="build"
    generator="Ninja"  # 可选 "Visual Studio 17 2022"（需手动修复路径）
    vs_tools_path_1=r"C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\Auxiliary\\Build\\vcvars64.bat"
    vs_tools_path_2=r"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat"

    if (os.path.exists(vs_tools_path_1)):
        vs_tools_path=vs_tools_path_1
    elif (os.path.exists(vs_tools_path_2)):
        vs_tools_path=vs_tools_path_2

    generate_compile_commands(
        project_root,
        build_dir,
        generator,
        vs_tools_path,
    )