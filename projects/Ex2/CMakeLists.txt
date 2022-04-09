add_executable(hf hf.c util.c)

target_compile_options(hf PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(hf m myfunc GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex2-hf
    COMMAND hf
    DEPENDS hf 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)