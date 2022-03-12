add_executable(bound bound.c util.c)

target_compile_options(bound PRIVATE -Wall -O3 -ffast-math)

target_link_libraries(bound m myfunc)

# add custom commands
add_custom_target(run-ex1-bound
    COMMAND bound
    DEPENDS bound 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)


###################
add_executable(scatt scatt.c util.c)

target_compile_options(scatt PRIVATE -Wall -O3 -ffast-math)

target_link_libraries(scatt m myfunc)

# add custom commands
add_custom_target(run-ex1-scatt
    COMMAND scatt
    DEPENDS scatt 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)