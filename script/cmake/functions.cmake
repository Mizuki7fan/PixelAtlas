#cmakelists中需要用到的通用函数

#复制项目必要的文件
function(copy_file proj_name file_name source_dir target_dir)
    set (full_source_path "${source_dir}/${file_name}")
    set (full_target_path "${target_dir}/${file_name}") 
    message(STATUS "文件源位置: ${full_source_path}")
    message(STATUS "文件目标位置: ${full_target_path}")
    if (NOT EXISTS ${full_source_path})
    message(FATAL_ERROR ${full_source_path}"文件不存在")
    endif()
    if (NOT EXISTS ${full_target_path})
    add_custom_command(TARGET ${proj_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${full_source_path}
                ${full_target_path}
            COMMENT "将${file_name}从${source_dir}复制到${target_dir}"
        )
    endif()
endfunction()