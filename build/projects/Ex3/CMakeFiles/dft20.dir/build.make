# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.20.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.20.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/teresa/Desktop/Repositories/ComputationalPhysics2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build

# Include any dependencies generated for this target.
include projects/Ex3/CMakeFiles/dft20.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include projects/Ex3/CMakeFiles/dft20.dir/compiler_depend.make

# Include the progress variables for this target.
include projects/Ex3/CMakeFiles/dft20.dir/progress.make

# Include the compile flags for this target's objects.
include projects/Ex3/CMakeFiles/dft20.dir/flags.make

projects/Ex3/CMakeFiles/dft20.dir/dft20.c.o: projects/Ex3/CMakeFiles/dft20.dir/flags.make
projects/Ex3/CMakeFiles/dft20.dir/dft20.c.o: ../projects/Ex3/dft20.c
projects/Ex3/CMakeFiles/dft20.dir/dft20.c.o: projects/Ex3/CMakeFiles/dft20.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object projects/Ex3/CMakeFiles/dft20.dir/dft20.c.o"
	cd /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3 && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT projects/Ex3/CMakeFiles/dft20.dir/dft20.c.o -MF CMakeFiles/dft20.dir/dft20.c.o.d -o CMakeFiles/dft20.dir/dft20.c.o -c /Users/teresa/Desktop/Repositories/ComputationalPhysics2/projects/Ex3/dft20.c

projects/Ex3/CMakeFiles/dft20.dir/dft20.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/dft20.dir/dft20.c.i"
	cd /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3 && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/teresa/Desktop/Repositories/ComputationalPhysics2/projects/Ex3/dft20.c > CMakeFiles/dft20.dir/dft20.c.i

projects/Ex3/CMakeFiles/dft20.dir/dft20.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/dft20.dir/dft20.c.s"
	cd /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3 && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/teresa/Desktop/Repositories/ComputationalPhysics2/projects/Ex3/dft20.c -o CMakeFiles/dft20.dir/dft20.c.s

# Object files for target dft20
dft20_OBJECTS = \
"CMakeFiles/dft20.dir/dft20.c.o"

# External object files for target dft20
dft20_EXTERNAL_OBJECTS =

projects/Ex3/dft20: projects/Ex3/CMakeFiles/dft20.dir/dft20.c.o
projects/Ex3/dft20: projects/Ex3/CMakeFiles/dft20.dir/build.make
projects/Ex3/dft20: myfunc/libmyfunc.a
projects/Ex3/dft20: /opt/homebrew/Cellar/gsl/2.6/lib/libgsl.dylib
projects/Ex3/dft20: /opt/homebrew/Cellar/gsl/2.6/lib/libgslcblas.dylib
projects/Ex3/dft20: /opt/homebrew/Cellar/gsl/2.6/lib/libgsl.dylib
projects/Ex3/dft20: /opt/homebrew/Cellar/gsl/2.6/lib/libgslcblas.dylib
projects/Ex3/dft20: projects/Ex3/CMakeFiles/dft20.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable dft20"
	cd /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dft20.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
projects/Ex3/CMakeFiles/dft20.dir/build: projects/Ex3/dft20
.PHONY : projects/Ex3/CMakeFiles/dft20.dir/build

projects/Ex3/CMakeFiles/dft20.dir/clean:
	cd /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3 && $(CMAKE_COMMAND) -P CMakeFiles/dft20.dir/cmake_clean.cmake
.PHONY : projects/Ex3/CMakeFiles/dft20.dir/clean

projects/Ex3/CMakeFiles/dft20.dir/depend:
	cd /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/teresa/Desktop/Repositories/ComputationalPhysics2 /Users/teresa/Desktop/Repositories/ComputationalPhysics2/projects/Ex3 /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3 /Users/teresa/Desktop/Repositories/ComputationalPhysics2/build/projects/Ex3/CMakeFiles/dft20.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : projects/Ex3/CMakeFiles/dft20.dir/depend
