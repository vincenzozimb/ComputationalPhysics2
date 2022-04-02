###########################################################
add_executable(scatt scatt.c util.c)

target_compile_options(scatt PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(scatt m myfunc GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex1-scatt
    COMMAND scatt && python plot_scatt.py
    DEPENDS scatt 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_scatt.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


###########################################################
add_executable(bound1d bound1d.c util.c)

target_compile_options(bound1d PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(bound1d m myfunc GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex1-bound1d
    COMMAND bound1d && gnuplot plot_psi1d.gp
    DEPENDS bound1d 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_psi1d.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

###########################################################
add_executable(bound3d bound3d.c util.c)

target_compile_options(bound3d PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(bound3d m myfunc GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex1-bound3d
    COMMAND bound3d && gnuplot plot_psi3d.gp
    DEPENDS bound3d 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_psi3d.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

###########################################################
add_executable(bessel_test bessel_test.c)

target_compile_options(bessel_test PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(bessel_test m myfunc GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex1-bessel_test
    COMMAND bessel_test && python plot_bessel.py
    DEPENDS bessel_test 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_bessel.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

