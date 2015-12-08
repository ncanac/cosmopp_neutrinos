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
find_library(COSMOPPLIB cosmopp ${COSMOPP_DIR}/build/lib)
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

# mn_scan
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(mn_scan mn_scan.cpp)
#target_link_libraries(mn_scan ${COSMOPPLIB})
#target_link_libraries(mn_scan ${MPI_CXX_LIBRARIES})
#target_link_libraries(mn_scan ${LAPACK_LIB_FLAGS})
#target_link_libraries(mn_scan ${CLASSLIB})
#target_link_libraries(mn_scan ${PLANCKLIB})
#target_link_libraries(mn_scan ${GSL_LIBRARIES})
#target_link_libraries(mn_scan ${WMAP9LIB})
#target_link_libraries(mn_scan ${CFITSIOLIB})
#target_link_libraries(mn_scan ${MULTINESTLIB})

# mn_planck
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(example_mn_planck example_mn_planck.cpp)
target_link_libraries(example_mn_planck ${COSMOPPLIB})
target_link_libraries(example_mn_planck ${MPI_CXX_LIBRARIES})
target_link_libraries(example_mn_planck ${LAPACK_LIB_FLAGS})
target_link_libraries(example_mn_planck ${CLASSLIB})
target_link_libraries(example_mn_planck ${PLANCKLIB})
#target_link_libraries(example_mn_planck ${GSL_LIBRARIES})
#target_link_libraries(example_mn_planck ${WMAP9LIB})
target_link_libraries(example_mn_planck ${CFITSIOLIB})
target_link_libraries(example_mn_planck ${MULTINESTLIB})
