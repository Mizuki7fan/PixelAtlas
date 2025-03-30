import os
import subprocess
import shutil
import json

def get_project_root(current_file: str) -> str:
    """根据脚本路径动态推导项目根目录（脚本位于 /pa-script/etc/）"""
    script_dir = os.path.dirname(os.path.abspath(current_file))  # G:/PixelAtlas/pa-script/etc
    project_root = os.path.abspath(os.path.join(script_dir, "../../"))  # 上溯两级到项目根目录
    return project_root

if __name__ == "__main__":
    project_root = get_project_root(__file__)
    build_dir="build"
    vs_tools_path_1=r"C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\Auxiliary\\Build\\vcvars64.bat"
    vs_tools_path_2=r"C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat"

    if (os.path.exists(vs_tools_path_1)):
        vs_tools_path=vs_tools_path_1
    elif (os.path.exists(vs_tools_path_2)):
        vs_tools_path=vs_tools_path_2

    # 清理并创建构建目录
    build_path = os.path.join(project_root, build_dir)
    if os.path.exists(build_path):
        shutil.rmtree(build_path)
    os.makedirs(build_path, exist_ok=True)

    clang_path = r"D:/LLVM/bin"
    env=os.environ.copy()
    env["CC"]=os.path.join(clang_path, "clang.exe")
    env["CXX"]=os.path.join(clang_path, "clang++.exe")

    cmake_cmd = [
        "cmake",
        "-G", "Ninja",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-B", build_path,
        "-S", project_root
    ]

    if os.name=="nt":
        subprocess.run(f'"{vs_tools_path}" && cmake {" ".join(cmake_cmd[1:])}', shell=True, cwd=build_path, env=env)

 # 删除除compile_commands.json外的所有文件
    for filename in os.listdir(build_path):
        if filename == "compile_commands.json":
            continue
        filepath = os.path.join(build_path, filename)
        if os.path.isfile(filepath) or os.path.islink(filepath):
            os.unlink(filepath)
        else:
            shutil.rmtree(filepath)

    # cmake_cmd = [
    #     "cmake",
    #     "-G", "Ninja",
    #     "-B", build_path,
    #     "-S", project_root
    # ]

    # if os.name=="nt":
    #     subprocess.run(f'"{vs_tools_path}" && cmake {" ".join(cmake_cmd[1:])}', shell=True, cwd=build_path, env=env)