###########################################################
add_executable(scatt scatt.c util.c)

target_compile_options(scatt PRIVATE -Wall -O3)

target_link_libraries(scatt m myfunc gsl)

# add custom commands
add_custom_target(run-ex1-scatt
    COMMAND scatt && python plot_scatt.py
    DEPENDS scatt 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_scatt.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


###########################################################
add_executable(bound bound.c util.c)

target_compile_options(bound PRIVATE -Wall -O3)

target_link_libraries(bound m myfunc gsl)

# add custom commands
add_custom_target(run-ex1-bound
    COMMAND bound
    DEPENDS bound 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)


###########################################################
add_executable(bessel_test bessel_test.c)

target_compile_options(bessel_test PRIVATE -Wall -O3)

target_link_libraries(bessel_test m myfunc gsl)

# add custom commands
add_custom_target(run-ex1-bessel_test
    COMMAND bessel_test && python plot_bessel.py
    DEPENDS bessel_test 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_bessel.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

