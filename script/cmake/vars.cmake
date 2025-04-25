#包含项目需要包含的宏变量
set(ACTIONLIST_FILE_NAME "pipeline.json")
add_definitions(-DACTIONLIST_FILE_NAME="${ACTIONLIST_FILE_NAME}")
set(RUN_DIR_NAME "z-run")
add_definitions(-DRUN_DIR_NAME="${RUN_DIR_NAME}")

set(LINSYSSOLVER "Eigen")