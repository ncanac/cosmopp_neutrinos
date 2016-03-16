cmake_minimum_required (VERSION 2.8.10)

project("test" CXX)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE)
endif (NOT CMAKE_BUILD_TYPE)

# Set directories
set(CLASS_DIR /Volumes/Data1/class_public-2.4.2)
#set(CLASS_DIR /Users/ncanac/Codes/class_public)
#set(COSMOPP_DIR /Volumes/Data1/COSMO++/cosmo_pp_private)
set(COSMOPP_DIR /Users/ncanac/Codes/cosmo_pp_private)
set(PLANCK_DIR /Volumes/Data1/Planck2015/plc-2.0)
set(PLANCK_DATA_DIR /Volumes/Data1/Planck2015/plc-2.0)
set(PLANCK_2015 ON)
#set(WMAP9_DIR /Volumes/Data1/WMAP9/wmap_likelihood_v5)
set(WMAP9_DIR /Users/ncanac/Codes/wmap_likelihood_v5)
set(PROJECT_DIR /Volumes/Data1/ncanac/cosmopp_neutrinos) 
set(CFITSIO_DIR /usr/local/cfitsio)
set(MULTINEST_DIR /usr/local)

# compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb")
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O2 -ffast-math")
	#set(CMAKE_CXX_FLAGS_RELEASE "-O2")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
	set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
else(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

# check for MPI
find_package(MPI REQUIRED)
add_definitions(-DCOSMO_MPI)
#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
if(NOT MPI_CXX_LIBRARIES)
	message(FATAL_ERROR "MPI found but the MPI C++ libraries not found!")
endif(NOT MPI_CXX_LIBRARIES)

# lapack
set(LAPACK_LIB_FLAGS "-framework Accelerate")
add_definitions(-DCOSMO_LAPACK)

# cfitsio
include_directories(${CFITSIO_DIR}/include)
find_library(CFITSIOLIB cfitsio ${CFITSIO_DIR}/lib)
if(NOT CFITSIOLIB)
	message(FATAL_ERROR "cfitsio library not found!")
endif(NOT CFITSIOLIB)

# planck
include_directories(${PLANCK_DIR}/include)
add_definitions(-DCOSMO_PLANCK -DPLANCK_DATA_DIR=${PLANCK_DATA_DIR})
if(PLANCK_2015)
	add_definitions(-DCOSMO_PLANCK_15)
endif(PLANCK_2015)
if(PLANCK_COMPILE_FLAGS)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PLANCK_COMPILE_FLAGS}")
endif(PLANCK_COMPILE_FLAGS)
find_library(PLANCKLIB clik ${PLANCK_DIR}/lib)
if(NOT PLANCKLIB)
	message(FATAL_ERROR "clik library not found!")
endif(NOT PLANCKLIB)

# WMAP
#include_directories(${WMAP9_DIR})
set(USE_FORT TRUE)
enable_language(Fortran)
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
add_definitions(-DWMAP9_IFORT)
add_definitions(-DCOSMO_WMAP9)
find_library(WMAP9LIB wmap9 ${WMAP9_DIR} NO_DEFAULT_PATH)
if(NOT WMAP9LIB)
    message(FATAL_ERROR "WMAP9 library not found!")
endif(NOT WMAP9LIB)

# class
include_directories(${CLASS_DIR}/include)
add_definitions(-DCOSMO_CLASS)

find_library(CLASSLIB class ${CLASS_DIR})
if(NOT CLASSLIB)
	message(FATAL_ERROR "class library not found!")
endif(NOT CLASSLIB)

#multinest
if(MULTINEST_DIR)
	include_directories(${MULTINEST_DIR}/include)
	add_definitions(-DCOSMO_MULTINEST)

	if(MPI_FOUND)
		find_library(MULTINESTLIB multinest_mpi ${MULTINEST_DIR}/lib)
	else(MPI_FOUND)
		find_library(MULTINESTLIB multinest ${MULTINEST_DIR}/lib)
	endif(MPI_FOUND)

	if(NOT MULTINESTLIB)
		message(FATAL_ERROR "multinest library not found!")
	endif(NOT MULTINESTLIB)
endif(MULTINEST_DIR)

# Cosmo++
find_library(COSMOPPLIB cosmopp ${COSMOPP_DIR}/lib)
include_directories(${COSMOPP_DIR}/include)

# openmp
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_definitions(-DCOSMO_OMP)

# GSL
find_package(GSL REQUIRED)
include_directories(GSL_INCLUDE_DIRS)

# Add project directory
include_directories(${PROJECT_DIR})

# Definitions
add_definitions(-DOUTPUT_LOG)
add_definitions(-DVERBOSE1)
add_definitions(-DCHECKS_ON)

# compile combined likelihood test 
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_combined_like test_combined_like.cpp sample.f90)
#target_link_libraries(test_combined_like ${COSMOPPLIB})
#target_link_libraries(test_combined_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_combined_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_combined_like ${CLASSLIB})
#target_link_libraries(test_combined_like ${PLANCKLIB})
#target_link_libraries(test_combined_like ${CLASSLIB})
#target_link_libraries(test_combined_like ${GSL_LIBRARIES})
#target_link_libraries(test_combined_like ${WMAP9LIB})
#target_link_libraries(test_combined_like ${CFITSIOLIB})

# compile learn_matrix
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(learn_matrix learn_matrix.cpp)
#target_link_libraries(learn_matrix ${COSMOPPLIB})
#target_link_libraries(learn_matrix ${MPI_CXX_LIBRARIES})
#target_link_libraries(learn_matrix ${LAPACK_LIB_FLAGS})

# compile BAO likelihood test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_bao_like test_bao_like.cpp)
#target_link_libraries(test_bao_like ${COSMOPPLIB})
#target_link_libraries(test_bao_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_bao_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_bao_like ${CLASSLIB})
#target_link_libraries(test_bao_like ${GSL_LIBRARIES})

# compile WiggleZ likelihood test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_wigglez_like test_wigglez_like.cpp)
#target_link_libraries(test_wigglez_like ${COSMOPPLIB})
#target_link_libraries(test_wigglez_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_wigglez_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_wigglez_like ${CLASSLIB})
#target_link_libraries(test_wigglez_like ${GSL_LIBRARIES})

# compile LRGHaloPs test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_LRGHaloPs test_LRGHaloPs.cpp)
#target_link_libraries(test_LRGHaloPs ${COSMOPPLIB})
#target_link_libraries(test_LRGHaloPs ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_LRGHaloPs ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_LRGHaloPs ${CLASSLIB})
#target_link_libraries(test_LRGHaloPs ${GSL_LIBRARIES})

# compile ps test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_ps test_ps.cpp)
#target_link_libraries(test_ps ${COSMOPPLIB})
#target_link_libraries(test_ps ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_ps ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_ps ${CLASSLIB})
#target_link_libraries(test_ps ${GSL_LIBRARIES})

# compile WMAP9 test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_wmap_like test_wmap_like.cpp sample.f90)
#target_link_libraries(test_wmap_like ${COSMOPPLIB})
#target_link_libraries(test_wmap_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_wmap_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_wmap_like ${CLASSLIB})
#target_link_libraries(test_wmap_like ${GSL_LIBRARIES})
#target_link_libraries(test_wmap_like ${WMAP9LIB})
#target_link_libraries(test_wmap_like ${CFITSIOLIB})

# compile sigma8 test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_sigma8 test_sigma8.cpp)
#target_link_libraries(test_sigma8 ${COSMOPPLIB})
#target_link_libraries(test_sigma8 ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_sigma8 ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_sigma8 ${CLASSLIB})

# compiles knotted multinest test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_knotted_planck_mn test_knotted_planck_mn.cpp sample.f90)
#target_link_libraries(test_knotted_planck_mn ${COSMOPPLIB})
#target_link_libraries(test_knotted_planck_mn ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_knotted_planck_mn ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_knotted_planck_mn ${CLASSLIB})
#target_link_libraries(test_knotted_planck_mn ${PLANCKLIB})
#target_link_libraries(test_knotted_planck_mn ${CFITSIOLIB})
#target_link_libraries(test_knotted_planck_mn ${MULTINESTLIB})

# compile LRG DR7 likelihood test
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(test_lrgdr7_like test_lrgdr7_like.cpp)
target_link_libraries(test_lrgdr7_like ${COSMOPPLIB})
target_link_libraries(test_lrgdr7_like ${MPI_CXX_LIBRARIES})
target_link_libraries(test_lrgdr7_like ${LAPACK_LIB_FLAGS})
target_link_libraries(test_lrgdr7_like ${CLASSLIB})
target_link_libraries(test_lrgdr7_like ${GSL_LIBRARIES})

# compiles knotted multinest test
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_knotted_nu_mn test_knotted_nu_mn.cpp sample.f90)
#target_link_libraries(test_knotted_nu_mn ${COSMOPPLIB})
#target_link_libraries(test_knotted_nu_mn ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_knotted_nu_mn ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_knotted_nu_mn ${CLASSLIB})
#target_link_libraries(test_knotted_nu_mn ${PLANCKLIB})
#target_link_libraries(test_knotted_nu_mn ${CFITSIOLIB})
#target_link_libraries(test_knotted_nu_mn ${MULTINESTLIB})
#target_link_libraries(test_knotted_nu_mn ${GSL_LIBRARIES})
#target_link_libraries(test_knotted_nu_mn ${WMAP9LIB})
