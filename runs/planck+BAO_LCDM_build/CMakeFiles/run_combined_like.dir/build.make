# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_SOURCE_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/runs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/planck+BAO_LCDM_build

# Include any dependencies generated for this target.
include CMakeFiles/run_combined_like.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/run_combined_like.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/run_combined_like.dir/flags.make

CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o: CMakeFiles/run_combined_like.dir/flags.make
CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o: ../run_combined_like.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/planck+BAO_LCDM_build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o"
	/sw/bin/g++-fsf-5   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o -c /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/run_combined_like.cpp

CMakeFiles/run_combined_like.dir/run_combined_like.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run_combined_like.dir/run_combined_like.cpp.i"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_FLAGS) -E /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/run_combined_like.cpp > CMakeFiles/run_combined_like.dir/run_combined_like.cpp.i

CMakeFiles/run_combined_like.dir/run_combined_like.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run_combined_like.dir/run_combined_like.cpp.s"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_FLAGS) -S /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/run_combined_like.cpp -o CMakeFiles/run_combined_like.dir/run_combined_like.cpp.s

CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.requires:
.PHONY : CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.requires

CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.provides: CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.requires
	$(MAKE) -f CMakeFiles/run_combined_like.dir/build.make CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.provides.build
.PHONY : CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.provides

CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.provides.build: CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o

# Object files for target run_combined_like
run_combined_like_OBJECTS = \
"CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o"

# External object files for target run_combined_like
run_combined_like_EXTERNAL_OBJECTS =

run_combined_like: CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o
run_combined_like: CMakeFiles/run_combined_like.dir/build.make
run_combined_like: /Volumes/Data1/COSMO++/cosmo_pp_private/build/lib/libcosmopp.a
run_combined_like: /usr/local/lib/libmpi_cxx.dylib
run_combined_like: /usr/local/lib/libmpi.dylib
run_combined_like: /usr/local/lib/libclass.a
run_combined_like: /usr/local/lib/libclik.dylib
run_combined_like: CMakeFiles/run_combined_like.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable run_combined_like"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run_combined_like.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/run_combined_like.dir/build: run_combined_like
.PHONY : CMakeFiles/run_combined_like.dir/build

CMakeFiles/run_combined_like.dir/requires: CMakeFiles/run_combined_like.dir/run_combined_like.cpp.o.requires
.PHONY : CMakeFiles/run_combined_like.dir/requires

CMakeFiles/run_combined_like.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run_combined_like.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run_combined_like.dir/clean

CMakeFiles/run_combined_like.dir/depend:
	cd /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/planck+BAO_LCDM_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Volumes/Data1/ncanac/cosmopp_neutrinos/runs /Volumes/Data1/ncanac/cosmopp_neutrinos/runs /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/planck+BAO_LCDM_build /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/planck+BAO_LCDM_build /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/planck+BAO_LCDM_build/CMakeFiles/run_combined_like.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/run_combined_like.dir/depend

