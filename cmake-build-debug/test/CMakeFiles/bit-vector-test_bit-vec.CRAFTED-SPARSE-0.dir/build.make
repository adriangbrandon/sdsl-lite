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

# Utility rule file for bit-vector-test_bit-vec.CRAFTED-SPARSE-0.

# Include the progress variables for this target.
include test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/progress.make

test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0: ../test/test_cases/bit-vec.CRAFTED-SPARSE-0
test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0: ../test/tmp
test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0: test/bit_vector_test
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Execute bit-vector-test on bit-vec.CRAFTED-SPARSE-0."
	/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/bit_vector_test /Users/adrian/CLionProjects/sdsl-lite/test/test_cases/bit-vec.CRAFTED-SPARSE-0 /Users/adrian/CLionProjects/sdsl-lite/test/tmp/bit_vector_test_bit-vec.CRAFTED-SPARSE-0 /Users/adrian/CLionProjects/sdsl-lite/test/tmp

../test/test_cases/bit-vec.CRAFTED-SPARSE-0: test/bit_vector_generator
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating test case bit-vec.CRAFTED-SPARSE-0."
	/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/bit_vector_generator /Users/adrian/CLionProjects/sdsl-lite/test/test_cases/bit-vec.CRAFTED-SPARSE-0 CRAFTED-SPARSE-0

bit-vector-test_bit-vec.CRAFTED-SPARSE-0: test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0
bit-vector-test_bit-vec.CRAFTED-SPARSE-0: ../test/test_cases/bit-vec.CRAFTED-SPARSE-0
bit-vector-test_bit-vec.CRAFTED-SPARSE-0: test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/build.make

.PHONY : bit-vector-test_bit-vec.CRAFTED-SPARSE-0

# Rule to build all files generated by this target.
test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/build: bit-vector-test_bit-vec.CRAFTED-SPARSE-0

.PHONY : test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/build

test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/clean:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && $(CMAKE_COMMAND) -P CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/clean

test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/depend:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/adrian/CLionProjects/sdsl-lite /Users/adrian/CLionProjects/sdsl-lite/test /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/bit-vector-test_bit-vec.CRAFTED-SPARSE-0.dir/depend

