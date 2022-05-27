add_executable(free free.c util.c)

target_compile_options(free PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(free m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex3
    COMMAND free && gnuplot plot.gp
    DEPENDS free 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})