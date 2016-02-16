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
CMAKE_SOURCE_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/runs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn

# Include any dependencies generated for this target.
include CMakeFiles/run_LCDM_mn.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/run_LCDM_mn.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/run_LCDM_mn.dir/flags.make

CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o: CMakeFiles/run_LCDM_mn.dir/flags.make
CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o: ../run_LCDM_mn.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o"
	/sw/bin/g++-fsf-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o -c /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/run_LCDM_mn.cpp

CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.i"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/run_LCDM_mn.cpp > CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.i

CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.s"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/run_LCDM_mn.cpp -o CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.s

CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.requires:

.PHONY : CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.requires

CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.provides: CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.requires
	$(MAKE) -f CMakeFiles/run_LCDM_mn.dir/build.make CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.provides.build
.PHONY : CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.provides

CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.provides.build: CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o


CMakeFiles/run_LCDM_mn.dir/sample.f90.o: CMakeFiles/run_LCDM_mn.dir/flags.make
CMakeFiles/run_LCDM_mn.dir/sample.f90.o: ../sample.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/run_LCDM_mn.dir/sample.f90.o"
	/usr/bin/ifort  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/sample.f90 -o CMakeFiles/run_LCDM_mn.dir/sample.f90.o

CMakeFiles/run_LCDM_mn.dir/sample.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/run_LCDM_mn.dir/sample.f90.i"
	/usr/bin/ifort  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/sample.f90 > CMakeFiles/run_LCDM_mn.dir/sample.f90.i

CMakeFiles/run_LCDM_mn.dir/sample.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/run_LCDM_mn.dir/sample.f90.s"
	/usr/bin/ifort  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/sample.f90 -o CMakeFiles/run_LCDM_mn.dir/sample.f90.s

CMakeFiles/run_LCDM_mn.dir/sample.f90.o.requires:

.PHONY : CMakeFiles/run_LCDM_mn.dir/sample.f90.o.requires

CMakeFiles/run_LCDM_mn.dir/sample.f90.o.provides: CMakeFiles/run_LCDM_mn.dir/sample.f90.o.requires
	$(MAKE) -f CMakeFiles/run_LCDM_mn.dir/build.make CMakeFiles/run_LCDM_mn.dir/sample.f90.o.provides.build
.PHONY : CMakeFiles/run_LCDM_mn.dir/sample.f90.o.provides

CMakeFiles/run_LCDM_mn.dir/sample.f90.o.provides.build: CMakeFiles/run_LCDM_mn.dir/sample.f90.o


# Object files for target run_LCDM_mn
run_LCDM_mn_OBJECTS = \
"CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o" \
"CMakeFiles/run_LCDM_mn.dir/sample.f90.o"

# External object files for target run_LCDM_mn
run_LCDM_mn_EXTERNAL_OBJECTS =

run_LCDM_mn: CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o
run_LCDM_mn: CMakeFiles/run_LCDM_mn.dir/sample.f90.o
run_LCDM_mn: CMakeFiles/run_LCDM_mn.dir/build.make
run_LCDM_mn: /Users/ncanac/Codes/cosmo_pp_private/lib/libcosmopp.a
run_LCDM_mn: /usr/local/lib/libmpi_cxx.dylib
run_LCDM_mn: /usr/local/lib/libmpi.dylib
run_LCDM_mn: /usr/local/lib/libclass.a
run_LCDM_mn: /usr/local/lib/libclik.dylib
run_LCDM_mn: /sw/lib/libgsl.dylib
run_LCDM_mn: /sw/lib/libgslcblas.dylib
run_LCDM_mn: /Users/ncanac/Codes/wmap_likelihood_v5/libwmap9.a
run_LCDM_mn: /usr/local/lib/libcfitsio.dylib
run_LCDM_mn: /usr/local/lib/libmultinest_mpi.dylib
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libifport.a
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libifcore.a
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libimf.a
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libsvml.a
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libipgo.a
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libirc.a
run_LCDM_mn: /usr/bin/ifort-15.0-base/compiler/lib/libsvml.a
run_LCDM_mn: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/7.0.2/lib/darwin/libclang_rt.osx.a
run_LCDM_mn: CMakeFiles/run_LCDM_mn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable run_LCDM_mn"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run_LCDM_mn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/run_LCDM_mn.dir/build: run_LCDM_mn

.PHONY : CMakeFiles/run_LCDM_mn.dir/build

CMakeFiles/run_LCDM_mn.dir/requires: CMakeFiles/run_LCDM_mn.dir/run_LCDM_mn.cpp.o.requires
CMakeFiles/run_LCDM_mn.dir/requires: CMakeFiles/run_LCDM_mn.dir/sample.f90.o.requires

.PHONY : CMakeFiles/run_LCDM_mn.dir/requires

CMakeFiles/run_LCDM_mn.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run_LCDM_mn.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run_LCDM_mn.dir/clean

CMakeFiles/run_LCDM_mn.dir/depend:
	cd /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Volumes/Data1/ncanac/cosmopp_neutrinos/runs /Volumes/Data1/ncanac/cosmopp_neutrinos/runs /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn /Volumes/Data1/ncanac/cosmopp_neutrinos/runs/LCDM_lrg_mn/CMakeFiles/run_LCDM_mn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/run_LCDM_mn.dir/depend

