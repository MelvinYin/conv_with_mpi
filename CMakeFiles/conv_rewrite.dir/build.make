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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/melvinyin/Desktop/work/conv_with_mpi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/melvinyin/Desktop/work/conv_with_mpi

# Include any dependencies generated for this target.
include CMakeFiles/conv_rewrite.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/conv_rewrite.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/conv_rewrite.dir/flags.make

CMakeFiles/conv_rewrite.dir/main.cpp.o: CMakeFiles/conv_rewrite.dir/flags.make
CMakeFiles/conv_rewrite.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/melvinyin/Desktop/work/conv_with_mpi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/conv_rewrite.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/conv_rewrite.dir/main.cpp.o -c /Users/melvinyin/Desktop/work/conv_with_mpi/main.cpp

CMakeFiles/conv_rewrite.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/conv_rewrite.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/melvinyin/Desktop/work/conv_with_mpi/main.cpp > CMakeFiles/conv_rewrite.dir/main.cpp.i

CMakeFiles/conv_rewrite.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/conv_rewrite.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/melvinyin/Desktop/work/conv_with_mpi/main.cpp -o CMakeFiles/conv_rewrite.dir/main.cpp.s

# Object files for target conv_rewrite
conv_rewrite_OBJECTS = \
"CMakeFiles/conv_rewrite.dir/main.cpp.o"

# External object files for target conv_rewrite
conv_rewrite_EXTERNAL_OBJECTS =

conv_rewrite: CMakeFiles/conv_rewrite.dir/main.cpp.o
conv_rewrite: CMakeFiles/conv_rewrite.dir/build.make
conv_rewrite: /usr/local/Cellar/open-mpi/4.0.3/lib/libmpi.dylib
conv_rewrite: CMakeFiles/conv_rewrite.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/melvinyin/Desktop/work/conv_with_mpi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable conv_rewrite"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/conv_rewrite.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/conv_rewrite.dir/build: conv_rewrite

.PHONY : CMakeFiles/conv_rewrite.dir/build

CMakeFiles/conv_rewrite.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/conv_rewrite.dir/cmake_clean.cmake
.PHONY : CMakeFiles/conv_rewrite.dir/clean

CMakeFiles/conv_rewrite.dir/depend:
	cd /Users/melvinyin/Desktop/work/conv_with_mpi && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/melvinyin/Desktop/work/conv_with_mpi /Users/melvinyin/Desktop/work/conv_with_mpi /Users/melvinyin/Desktop/work/conv_with_mpi /Users/melvinyin/Desktop/work/conv_with_mpi /Users/melvinyin/Desktop/work/conv_with_mpi/CMakeFiles/conv_rewrite.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/conv_rewrite.dir/depend

