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
CMAKE_BINARY_DIR = /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build

# Include any dependencies generated for this target.
include CMakeFiles/test_combined_like.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_combined_like.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_combined_like.dir/flags.make

CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o: CMakeFiles/test_combined_like.dir/flags.make
CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o: ../test_combined_like.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o"
	/sw/bin/g++-fsf-5   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o -c /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/test_combined_like.cpp

CMakeFiles/test_combined_like.dir/test_combined_like.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_combined_like.dir/test_combined_like.cpp.i"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/test_combined_like.cpp > CMakeFiles/test_combined_like.dir/test_combined_like.cpp.i

CMakeFiles/test_combined_like.dir/test_combined_like.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_combined_like.dir/test_combined_like.cpp.s"
	/sw/bin/g++-fsf-5  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/test_combined_like.cpp -o CMakeFiles/test_combined_like.dir/test_combined_like.cpp.s

CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.requires:

.PHONY : CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.requires

CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.provides: CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_combined_like.dir/build.make CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.provides.build
.PHONY : CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.provides

CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.provides.build: CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o


CMakeFiles/test_combined_like.dir/sample.f90.o: CMakeFiles/test_combined_like.dir/flags.make
CMakeFiles/test_combined_like.dir/sample.f90.o: ../sample.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/test_combined_like.dir/sample.f90.o"
	/usr/bin/ifort  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sample.f90 -o CMakeFiles/test_combined_like.dir/sample.f90.o

CMakeFiles/test_combined_like.dir/sample.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/test_combined_like.dir/sample.f90.i"
	/usr/bin/ifort  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sample.f90 > CMakeFiles/test_combined_like.dir/sample.f90.i

CMakeFiles/test_combined_like.dir/sample.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/test_combined_like.dir/sample.f90.s"
	/usr/bin/ifort  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/sample.f90 -o CMakeFiles/test_combined_like.dir/sample.f90.s

CMakeFiles/test_combined_like.dir/sample.f90.o.requires:

.PHONY : CMakeFiles/test_combined_like.dir/sample.f90.o.requires

CMakeFiles/test_combined_like.dir/sample.f90.o.provides: CMakeFiles/test_combined_like.dir/sample.f90.o.requires
	$(MAKE) -f CMakeFiles/test_combined_like.dir/build.make CMakeFiles/test_combined_like.dir/sample.f90.o.provides.build
.PHONY : CMakeFiles/test_combined_like.dir/sample.f90.o.provides

CMakeFiles/test_combined_like.dir/sample.f90.o.provides.build: CMakeFiles/test_combined_like.dir/sample.f90.o


# Object files for target test_combined_like
test_combined_like_OBJECTS = \
"CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o" \
"CMakeFiles/test_combined_like.dir/sample.f90.o"

# External object files for target test_combined_like
test_combined_like_EXTERNAL_OBJECTS =

test_combined_like: CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o
test_combined_like: CMakeFiles/test_combined_like.dir/sample.f90.o
test_combined_like: CMakeFiles/test_combined_like.dir/build.make
test_combined_like: /Users/ncanac/Codes/cosmo_pp_private/lib/libcosmopp.a
test_combined_like: /usr/local/lib/libmpi_cxx.dylib
test_combined_like: /usr/local/lib/libmpi.dylib
test_combined_like: /usr/local/lib/libclass.a
test_combined_like: /usr/local/lib/libclik.dylib
test_combined_like: /usr/local/lib/libclass.a
test_combined_like: /sw/lib/libgsl.dylib
test_combined_like: /sw/lib/libgslcblas.dylib
test_combined_like: /Users/ncanac/Codes/wmap_likelihood_v5/libwmap9.a
test_combined_like: /usr/local/lib/libcfitsio.dylib
test_combined_like: /usr/local/lib/libclik.dylib
test_combined_like: /sw/lib/libgsl.dylib
test_combined_like: /sw/lib/libgslcblas.dylib
test_combined_like: /Users/ncanac/Codes/wmap_likelihood_v5/libwmap9.a
test_combined_like: /usr/local/lib/libcfitsio.dylib
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libifport.a
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libifcore.a
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libimf.a
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libsvml.a
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libipgo.a
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libirc.a
test_combined_like: /usr/bin/ifort-15.0-base/compiler/lib/libsvml.a
test_combined_like: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/7.0.2/lib/darwin/libclang_rt.osx.a
test_combined_like: CMakeFiles/test_combined_like.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable test_combined_like"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_combined_like.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_combined_like.dir/build: test_combined_like

.PHONY : CMakeFiles/test_combined_like.dir/build

CMakeFiles/test_combined_like.dir/requires: CMakeFiles/test_combined_like.dir/test_combined_like.cpp.o.requires
CMakeFiles/test_combined_like.dir/requires: CMakeFiles/test_combined_like.dir/sample.f90.o.requires

.PHONY : CMakeFiles/test_combined_like.dir/requires

CMakeFiles/test_combined_like.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_combined_like.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_combined_like.dir/clean

CMakeFiles/test_combined_like.dir/depend:
	cd /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Volumes/Data1/ncanac/cosmopp_neutrinos/tests /Volumes/Data1/ncanac/cosmopp_neutrinos/tests /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build /Volumes/Data1/ncanac/cosmopp_neutrinos/tests/build/CMakeFiles/test_combined_like.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_combined_like.dir/depend
