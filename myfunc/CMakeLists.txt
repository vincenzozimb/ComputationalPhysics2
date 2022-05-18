# Adding the library
add_library(myfunc 
                    print_routines.c
                    numerov.c
                    bisection.c
                    zeros_newton.c
                    integrate_trap.c
                    bessel_func.c
                    array_routines.c
                    lapack_wrappers.c
)

target_compile_options(myfunc PRIVATE -Wall -O3)
target_link_libraries(myfunc m)