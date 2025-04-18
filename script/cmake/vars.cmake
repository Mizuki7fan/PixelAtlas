#包含项目需要包含的宏变量
set(ACTIONLIST_FILE "pipeline.json")
add_definitions(-DACTIONLIST_FILE="${ACTIONLIST_FILE}")

set(LINSYSSOLVER "Eigen")