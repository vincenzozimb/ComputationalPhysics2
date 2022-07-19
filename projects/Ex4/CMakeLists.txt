###########################################################
add_executable(mc1d mc1d.c )

target_compile_options(mc1d PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(mc1d m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex4-mc1d
COMMAND mc1d && gnuplot plot.gp
DEPENDS mc1d 
WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

###########################################################
add_executable(mc3d mc3d.c )

target_compile_options(mc3d PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(mc3d m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex4-mc3d
COMMAND mc3d && gnuplot plot3d.gp
DEPENDS mc3d 
WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot3d.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

###########################################################
add_executable(mcel mcel.c util.c)

target_compile_options(mcel PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(mcel m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex4-mcel
COMMAND mcel && gnuplot plotEl.gp
DEPENDS mcel 
WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plotEl.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})



