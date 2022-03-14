add_executable(test test.c util.c)

target_compile_options(test PRIVATE -Wall -O3)

target_link_libraries(test m myfunc)

# add custom commands
add_custom_target(run-test
    COMMAND test && python plot.py
    DEPENDS test 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})