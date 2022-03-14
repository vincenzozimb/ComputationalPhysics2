add_executable(test test.c)

target_compile_options(test PRIVATE -Wall -O3)

target_link_libraries(test m myfunc)

# add custom commands
add_custom_target(run-test
    COMMAND test
    DEPENDS test 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)