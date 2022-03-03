add_library(myCP-lib STATIC
                        numerical_methods/integrate.c
                        numerical_methods/minimum.c
                        numerical_methods/zeros_newton.c
                        numerical_methods/bisection.c
                        numerical_methods/derivate.c
                        #numerical_methods/fast_squareroot.c
                        differential_eq/numerov.c
                        differential_eq/verlet.c
                        differential_eq/runge_kutta.c
                        linear_algebra/array_routines.c
                        linear_algebra/blas_wrappers.c
                        linear_algebra/lapack_wrappers.c
                        FFT_routines/radial.c
                        FFT_routines/direct_wrappers.c
                        print_routines.c
                        gnuplot_i.c
)

target_include_directories(myCP-lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_options(myCP-lib PRIVATE -Wall -O3 -march=native -ffast-math -std=c11)
target_link_libraries(myCP-lib)

#target_link_libraries(myCP-lib extern/gnuplot_i)

# Collect header files in include dir for installation
# file(GLOB_RECURSE INCLUDE_FILES "${PROJECT_SOURCE_DIR}/include/*.h")
# set_target_properties(myCP-lib PROPERTIES PUBLIC_HEADER "${INCLUDE_FILES}")

#include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR $ENV{HOME}/.local)
install( 
    TARGETS myCP-lib
    EXPORT MyCPLibTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_LIBDIR}/include
)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    FILES_MATCHING PATTERN "*.h"
)

install(
    EXPORT MyCPLibTargets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MyCPLib
)

include(CMakePackageConfigHelpers)
configure_package_config_file( 
    "Config.cmake.in" 
    "MyCPLibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MyCPLib
    PATH_VARS
    CMAKE_INSTALL_LIBDIR
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MyCPLibConfigVersion.cmake
    VERSION 1.0.0
    COMPATIBILITY SameMajorVersion
)

### Install Config and ConfigVersion files
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/MyCPLibConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/MyCPLibConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/MyCPLib"
)
