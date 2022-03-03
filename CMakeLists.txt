cmake_minimum_required(VERSION 3.13)
project(ComputationalPhysics2)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)


# List here the directories in extern
add_subdirectory(extern/MyCP_libraries)


# List here the directories in src
add_subdirectory(src/test)
