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
include external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/depend.make

# Include the progress variables for this target.
include external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/progress.make

# Include the compile flags for this target's objects.
include external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/flags.make

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/divsufsort.o: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/flags.make
external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/divsufsort.o: ../external/libdivsufsort/lib/divsufsort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/divsufsort.o"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort.dir/divsufsort.o   -c /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/divsufsort.c

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/divsufsort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort.dir/divsufsort.i"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/divsufsort.c > CMakeFiles/divsufsort.dir/divsufsort.i

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/divsufsort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort.dir/divsufsort.s"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/divsufsort.c -o CMakeFiles/divsufsort.dir/divsufsort.s

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/sssort.o: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/flags.make
external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/sssort.o: ../external/libdivsufsort/lib/sssort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/sssort.o"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort.dir/sssort.o   -c /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/sssort.c

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/sssort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort.dir/sssort.i"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/sssort.c > CMakeFiles/divsufsort.dir/sssort.i

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/sssort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort.dir/sssort.s"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/sssort.c -o CMakeFiles/divsufsort.dir/sssort.s

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/trsort.o: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/flags.make
external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/trsort.o: ../external/libdivsufsort/lib/trsort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/trsort.o"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort.dir/trsort.o   -c /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/trsort.c

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/trsort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort.dir/trsort.i"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/trsort.c > CMakeFiles/divsufsort.dir/trsort.i

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/trsort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort.dir/trsort.s"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/trsort.c -o CMakeFiles/divsufsort.dir/trsort.s

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/utils.o: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/flags.make
external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/utils.o: ../external/libdivsufsort/lib/utils.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/utils.o"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort.dir/utils.o   -c /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/utils.c

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/utils.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort.dir/utils.i"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/utils.c > CMakeFiles/divsufsort.dir/utils.i

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/utils.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort.dir/utils.s"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib/utils.c -o CMakeFiles/divsufsort.dir/utils.s

# Object files for target divsufsort
divsufsort_OBJECTS = \
"CMakeFiles/divsufsort.dir/divsufsort.o" \
"CMakeFiles/divsufsort.dir/sssort.o" \
"CMakeFiles/divsufsort.dir/trsort.o" \
"CMakeFiles/divsufsort.dir/utils.o"

# External object files for target divsufsort
divsufsort_EXTERNAL_OBJECTS =

external/libdivsufsort/lib/libdivsufsort.a: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/divsufsort.o
external/libdivsufsort/lib/libdivsufsort.a: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/sssort.o
external/libdivsufsort/lib/libdivsufsort.a: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/trsort.o
external/libdivsufsort/lib/libdivsufsort.a: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/utils.o
external/libdivsufsort/lib/libdivsufsort.a: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/build.make
external/libdivsufsort/lib/libdivsufsort.a: external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C static library libdivsufsort.a"
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && $(CMAKE_COMMAND) -P CMakeFiles/divsufsort.dir/cmake_clean_target.cmake
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/divsufsort.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/build: external/libdivsufsort/lib/libdivsufsort.a

.PHONY : external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/build

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/clean:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib && $(CMAKE_COMMAND) -P CMakeFiles/divsufsort.dir/cmake_clean.cmake
.PHONY : external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/clean

external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/depend:
	cd /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/adrian/CLionProjects/sdsl-lite /Users/adrian/CLionProjects/sdsl-lite/external/libdivsufsort/lib /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib /Users/adrian/CLionProjects/sdsl-lite/cmake-build-debug/external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/libdivsufsort/lib/CMakeFiles/divsufsort.dir/depend

