cmake_minimum_required (VERSION 2.8.10)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE)
endif (NOT CMAKE_BUILD_TYPE)

# Set directories
#set(COSMOPP_DIR /Volumes/Data1/COSMO++/cosmo_pp_private)
set(COSMOPP_DIR /Users/ncanac/Codes/cosmo_pp_private)
set(PROJECT_DIR /Volumes/Data1/ncanac/cosmopp_neutrinos) 

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


# Cosmo++
find_library(COSMOPPLIB cosmopp ${COSMOPP_DIR}/build/lib)
include_directories(${COSMOPP_DIR}/include)

# openmp
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_definitions(-DCOSMO_OMP)

# check for MPI
find_package(MPI REQUIRED)
add_definitions(-DCOSMO_MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
if(NOT MPI_CXX_LIBRARIES)
	message(FATAL_ERROR "MPI found but the MPI C++ libraries not found!")
endif(NOT MPI_CXX_LIBRARIES)

# lapack
set(LAPACK_LIB_FLAGS "-framework Accelerate")
add_definitions(-DCOSMO_LAPACK)

# Add project directory
include_directories(${PROJECT_DIR})

# Definitions
add_definitions(-DOUTPUT_LOG)
add_definitions(-DVERBOSE1)

# planck+bao 
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(analyze_MC analyze_MC.cpp)
target_link_libraries(analyze_MC ${COSMOPPLIB})
target_link_libraries(analyze_MC ${MPI_CXX_LIBRARIES})
target_link_libraries(analyze_MC ${LAPACK_LIB_FLAGS})
