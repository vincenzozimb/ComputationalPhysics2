###########################################################
add_executable(free free.c util.c)

target_compile_options(free PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(scatt m myfunc GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex3
    COMMAND free
    DEPENDS free 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)