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

# Include any dependencies generated for this target.
include test/CMakeFiles/wt_int_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/wt_int_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/wt_int_test.dir/flags.make

test/CMakeFiles/wt_int_test.dir/wt_int_test.cpp.o: test/CMakeFiles/wt_int_test.dir/flags.make
test/CMakeFiles/wt_int_test.dir/wt_int_test.cpp.o: ../test/wt_int_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/wt_int_test.dir/wt_int_test.cpp.o"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wt_int_test.dir/wt_int_test.cpp.o -c /Users/adrian/CLionProjects/sdsl-lite/test/wt_int_test.cpp

test/CMakeFiles/wt_int_test.dir/wt_int_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wt_int_test.dir/wt_int_test.cpp.i"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/adrian/CLionProjects/sdsl-lite/test/wt_int_test.cpp > CMakeFiles/wt_int_test.dir/wt_int_test.cpp.i

test/CMakeFiles/wt_int_test.dir/wt_int_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wt_int_test.dir/wt_int_test.cpp.s"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/adrian/CLionProjects/sdsl-lite/test/wt_int_test.cpp -o CMakeFiles/wt_int_test.dir/wt_int_test.cpp.s

# Object files for target wt_int_test
wt_int_test_OBJECTS = \
"CMakeFiles/wt_int_test.dir/wt_int_test.cpp.o"

# External object files for target wt_int_test
wt_int_test_EXTERNAL_OBJECTS =

test/wt_int_test: test/CMakeFiles/wt_int_test.dir/wt_int_test.cpp.o
test/wt_int_test: test/CMakeFiles/wt_int_test.dir/build.make
test/wt_int_test: lib/libsdsl.a
test/wt_int_test: external/googletest/googletest/libgtest.a
test/wt_int_test: external/libdivsufsort/lib/libdivsufsort.a
test/wt_int_test: external/libdivsufsort/lib/libdivsufsort64.a
test/wt_int_test: test/CMakeFiles/wt_int_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable wt_int_test"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wt_int_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/wt_int_test.dir/build: test/wt_int_test

.PHONY : test/CMakeFiles/wt_int_test.dir/build

test/CMakeFiles/wt_int_test.dir/clean:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test && $(CMAKE_COMMAND) -P CMakeFiles/wt_int_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/wt_int_test.dir/clean

test/CMakeFiles/wt_int_test.dir/depend:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/adrian/CLionProjects/sdsl-lite /Users/adrian/CLionProjects/sdsl-lite/test /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/test/CMakeFiles/wt_int_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/wt_int_test.dir/depend

