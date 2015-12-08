cmake_minimum_required (VERSION 2.8.10)

project("runs" CXX)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE)
endif (NOT CMAKE_BUILD_TYPE)

# Set directories
set(CLASS_DIR /Volumes/Data1/class_public-2.4.2)
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
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
	set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
else(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

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

# Add project directory
include_directories(${PROJECT_DIR})

# Definitions
add_definitions(-DOUTPUT_LOG)
add_definitions(-DVERBOSE1)

# bao
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_bao_like run_bao_like.cpp)
#target_link_libraries(run_bao_like ${COSMOPPLIB})
#target_link_libraries(run_bao_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_bao_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_bao_like ${CLASSLIB})

# planck
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_planck_like run_planck_like.cpp)
#target_link_libraries(run_planck_like ${COSMOPPLIB})
#target_link_libraries(run_planck_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_planck_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_planck_like ${CLASSLIB})
#target_link_libraries(run_planck_like ${PLANCKLIB})

# combined
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_combined_like run_combined_like.cpp)
#target_link_libraries(run_combined_like ${COSMOPPLIB})
#target_link_libraries(run_combined_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_combined_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_combined_like ${CLASSLIB})
#target_link_libraries(run_combined_like ${PLANCKLIB})

# LRG
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_lrg_like run_lrg_like.cpp)
#target_link_libraries(run_lrg_like ${COSMOPPLIB})
#target_link_libraries(run_lrg_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_lrg_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_lrg_like ${CLASSLIB})
#target_link_libraries(run_lrg_like ${GSL_LIBRARIES})

# WMAP + LRG
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_wmap+lrg_like run_wmap+lrg_like.cpp sample.f90)
#target_link_libraries(run_wmap+lrg_like ${COSMOPPLIB})
#target_link_libraries(run_wmap+lrg_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_wmap+lrg_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_wmap+lrg_like ${CLASSLIB})
#target_link_libraries(run_wmap+lrg_like ${GSL_LIBRARIES})
#target_link_libraries(run_wmap+lrg_like ${WMAP9LIB})
#target_link_libraries(run_wmap+lrg_like ${CFITSIOLIB})
#target_link_libraries(run_wmap+lrg_like ${PLANCKLIB})

# planck+bao 
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_planck+BAO_like run_planck+BAO_like.cpp sample.f90)
#target_link_libraries(run_planck+BAO_like ${COSMOPPLIB})
#target_link_libraries(run_planck+BAO_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_planck+BAO_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_planck+BAO_like ${CLASSLIB})
#target_link_libraries(run_planck+BAO_like ${PLANCKLIB})
#target_link_libraries(run_planck+BAO_like ${GSL_LIBRARIES})
#target_link_libraries(run_planck+BAO_like ${WMAP9LIB})
#target_link_libraries(run_planck+BAO_like ${CFITSIOLIB})

# knotted
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_knotted run_knotted.cpp sample.f90)
#target_link_libraries(run_knotted ${COSMOPPLIB})
#target_link_libraries(run_knotted ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_knotted ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_knotted ${CLASSLIB})
#target_link_libraries(run_knotted ${PLANCKLIB})
#target_link_libraries(run_knotted ${GSL_LIBRARIES})
#target_link_libraries(run_knotted ${WMAP9LIB})
#target_link_libraries(run_knotted ${CFITSIOLIB})
#target_link_libraries(run_knotted ${MULTINESTLIB})

# planck+lrg multinest
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_planck+lrg_mn_like run_planck+lrg_mn_like.cpp sample.f90)
#target_link_libraries(run_planck+lrg_mn_like ${COSMOPPLIB})
#target_link_libraries(run_planck+lrg_mn_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_planck+lrg_mn_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_planck+lrg_mn_like ${CLASSLIB})
#target_link_libraries(run_planck+lrg_mn_like ${PLANCKLIB})
#target_link_libraries(run_planck+lrg_mn_like ${GSL_LIBRARIES})
#target_link_libraries(run_planck+lrg_mn_like ${WMAP9LIB})
#target_link_libraries(run_planck+lrg_mn_like ${CFITSIOLIB})
#target_link_libraries(run_planck+lrg_mn_like ${MULTINESTLIB})

# knotted MCMC
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(run_knotted_MCMC run_knotted_MCMC.cpp sample.f90)
target_link_libraries(run_knotted_MCMC ${COSMOPPLIB})
target_link_libraries(run_knotted_MCMC ${MPI_CXX_LIBRARIES})
target_link_libraries(run_knotted_MCMC ${LAPACK_LIB_FLAGS})
target_link_libraries(run_knotted_MCMC ${CLASSLIB})
target_link_libraries(run_knotted_MCMC ${PLANCKLIB})
target_link_libraries(run_knotted_MCMC ${GSL_LIBRARIES})
target_link_libraries(run_knotted_MCMC ${WMAP9LIB})
target_link_libraries(run_knotted_MCMC ${CFITSIOLIB})
target_link_libraries(run_knotted_MCMC ${MULTINESTLIB})

# LCDM multinest
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_LCDM_mn run_LCDM_mn.cpp sample.f90)
#target_link_libraries(run_LCDM_mn ${COSMOPPLIB})
#target_link_libraries(run_LCDM_mn ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_LCDM_mn ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_LCDM_mn ${CLASSLIB})
#target_link_libraries(run_LCDM_mn ${PLANCKLIB})
#target_link_libraries(run_LCDM_mn ${GSL_LIBRARIES})
#target_link_libraries(run_LCDM_mn ${WMAP9LIB})
#target_link_libraries(run_LCDM_mn ${CFITSIOLIB})
#target_link_libraries(run_LCDM_mn ${MULTINESTLIB})

# standard PPS multinest
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_standard_mn run_standard_mn.cpp sample.f90)
#target_link_libraries(run_standard_mn ${COSMOPPLIB})
#target_link_libraries(run_standard_mn ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_standard_mn ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_standard_mn ${CLASSLIB})
#target_link_libraries(run_standard_mn ${PLANCKLIB})
#target_link_libraries(run_standard_mn ${GSL_LIBRARIES})
#target_link_libraries(run_standard_mn ${WMAP9LIB})
#target_link_libraries(run_standard_mn ${CFITSIOLIB})
#target_link_libraries(run_standard_mn ${MULTINESTLIB})

# knotted multinest no neutrinos
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(run_knotted_mn run_knotted_mn.cpp sample.f90)
#target_link_libraries(run_knotted_mn ${COSMOPPLIB})
#target_link_libraries(run_knotted_mn ${MPI_CXX_LIBRARIES})
#target_link_libraries(run_knotted_mn ${LAPACK_LIB_FLAGS})
#target_link_libraries(run_knotted_mn ${CLASSLIB})
#target_link_libraries(run_knotted_mn ${PLANCKLIB})
#target_link_libraries(run_knotted_mn ${GSL_LIBRARIES})
#target_link_libraries(run_knotted_mn ${WMAP9LIB})
#target_link_libraries(run_knotted_mn ${CFITSIOLIB})
#target_link_libraries(run_knotted_mn ${MULTINESTLIB})
