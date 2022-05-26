add_executable(dft dft.c)

target_compile_options(dft PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(dft m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex3
    COMMAND dft && gnuplot plot.gp
    DEPENDS dft 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})