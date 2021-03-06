# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.4

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /sw/bin/cmake

# The command to remove a file.
RM = /sw/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/tests

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test

# Include any dependencies generated for this target.
include CMakeFiles/test_sigma8.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_sigma8.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_sigma8.dir/flags.make

CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o: CMakeFiles/test_sigma8.dir/flags.make
CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o: ../test_sigma8.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o"
	/sw/bin/g++-fsf-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o -c /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/test_sigma8.cpp

CMakeFiles/test_sigma8.dir/test_sigma8.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_sigma8.dir/test_sigma8.cpp.i"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/test_sigma8.cpp > CMakeFiles/test_sigma8.dir/test_sigma8.cpp.i

CMakeFiles/test_sigma8.dir/test_sigma8.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_sigma8.dir/test_sigma8.cpp.s"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/test_sigma8.cpp -o CMakeFiles/test_sigma8.dir/test_sigma8.cpp.s

CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.requires:

.PHONY : CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.requires

CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.provides: CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_sigma8.dir/build.make CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.provides.build
.PHONY : CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.provides

CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.provides.build: CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o


# Object files for target test_sigma8
test_sigma8_OBJECTS = \
"CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o"

# External object files for target test_sigma8
test_sigma8_EXTERNAL_OBJECTS =

test_sigma8: CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o
test_sigma8: CMakeFiles/test_sigma8.dir/build.make
test_sigma8: /Users/ncanac/Codes/cosmo_pp_private/lib/libcosmopp.a
test_sigma8: /usr/local/lib/libmpi_cxx.dylib
test_sigma8: /usr/local/lib/libmpi.dylib
test_sigma8: /usr/local/lib/libclass.a
test_sigma8: CMakeFiles/test_sigma8.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_sigma8"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_sigma8.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_sigma8.dir/build: test_sigma8

.PHONY : CMakeFiles/test_sigma8.dir/build

CMakeFiles/test_sigma8.dir/requires: CMakeFiles/test_sigma8.dir/test_sigma8.cpp.o.requires

.PHONY : CMakeFiles/test_sigma8.dir/requires

CMakeFiles/test_sigma8.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_sigma8.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_sigma8.dir/clean

CMakeFiles/test_sigma8.dir/depend:
	cd /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Volumes/Data1/ncanac/cosmopp_neutrinos/tests /Volumes/Data1/ncanac/cosmopp_neutrinos/tests /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sigma8_test/CMakeFiles/test_sigma8.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_sigma8.dir/depend

