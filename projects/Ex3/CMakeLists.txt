########################################################### dft8
add_executable(dft8 dft8.c)

target_compile_options(dft8 PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(dft8 m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex3-8
COMMAND dft8 && gnuplot plot_wf_free.gp
DEPENDS dft8 
WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_wf_free.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

########################################################### dft20
add_executable(dft20 dft20.c)

target_compile_options(dft20 PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(dft20 m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex3-20
    COMMAND dft20 #&& gnuplot plot_20.gp
    DEPENDS dft20 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_20.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

########################################################### dft40
add_executable(dft40 dft40.c)

target_compile_options(dft40 PRIVATE -Wall -O3)
find_package(GSL REQUIRED)
target_link_libraries(dft40 m myfunc lapack GSL::gsl GSL::gslcblas)

# add custom commands
add_custom_target(run-ex3-40
    COMMAND dft40 #&& gnuplot plot_40.gp
    DEPENDS dft40 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
    )
    
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_40.gp
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    
########################################################### plot densities
# add custom commands
add_custom_target(plot-ex3-dft
    COMMAND gnuplot plot.gp
    DEPENDS plot.gp
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot.gp
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})