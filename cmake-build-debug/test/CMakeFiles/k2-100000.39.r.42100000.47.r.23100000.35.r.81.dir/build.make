# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/adrian/CLionProjects/sdsl-lite

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug

# Utility rule file for k2-100000.39.r.42100000.47.r.23100000.35.r.81.

# Include the progress variables for this target.
include test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/progress.make

test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81: ../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.x
test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81: ../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.y
test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81: ../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.w


../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.x: test/int_vector_generator
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating test case k2-100000.39.r.42100000.47.r.23100000.35.r.81."
	/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/int_vector_generator /Users/adrian/CLionProjects/sdsl-lite/test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.x 100000 39 r 42

../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.y: test/int_vector_generator
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating test case k2-100000.39.r.42100000.47.r.23100000.35.r.81."
	/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/int_vector_generator /Users/adrian/CLionProjects/sdsl-lite/test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.y 100000 47 r 23

../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.w: test/int_vector_generator
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating test case k2-100000.39.r.42100000.47.r.23100000.35.r.81."
	/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/int_vector_generator /Users/adrian/CLionProjects/sdsl-lite/test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.w 100000 35 r 81

k2-100000.39.r.42100000.47.r.23100000.35.r.81: test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81
k2-100000.39.r.42100000.47.r.23100000.35.r.81: ../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.x
k2-100000.39.r.42100000.47.r.23100000.35.r.81: ../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.y
k2-100000.39.r.42100000.47.r.23100000.35.r.81: ../test/test_cases/k2-100000.39.r.42100000.47.r.23100000.35.r.81.w
k2-100000.39.r.42100000.47.r.23100000.35.r.81: test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/build.make

.PHONY : k2-100000.39.r.42100000.47.r.23100000.35.r.81

# Rule to build all files generated by this target.
test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/build: k2-100000.39.r.42100000.47.r.23100000.35.r.81

.PHONY : test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/build

test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/clean:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && $(CMAKE_COMMAND) -P CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/clean

test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/depend:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/adrian/CLionProjects/sdsl-lite /Users/adrian/CLionProjects/sdsl-lite/test /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/k2-100000.39.r.42100000.47.r.23100000.35.r.81.dir/depend

