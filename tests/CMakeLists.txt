cmake_minimum_required (VERSION 2.8.10)

project("bao" CXX)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE)
endif (NOT CMAKE_BUILD_TYPE)

# Set directories
set(CLASS_DIR /Volumes/Data1/class_public-2.4.2)
set(COSMOPP_DIR /Volumes/Data1/COSMO++/cosmo_pp_private)
set(PLANCK_DIR /Volumes/Data1/Planck2015/plc-2.0)
set(PLANCK_DATA_DIR "/Volumes/Data1/Planck2015/plc-2.0")
set(PLANCK_2015 ON)
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

# class
include_directories(${CLASS_DIR}/include)
add_definitions(-DCOSMO_CLASS)

find_library(CLASSLIB class ${CLASS_DIR})
if(NOT CLASSLIB)
	message(FATAL_ERROR "class library not found!")
endif(NOT CLASSLIB)

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
#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
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
add_definitions(-DCHECKS_ON)

# compile combined likelihood test 
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_executable(test_combined_like test_combined_like.cpp)
#target_link_libraries(test_combined_like ${COSMOPPLIB})
#target_link_libraries(test_combined_like ${MPI_CXX_LIBRARIES})
#target_link_libraries(test_combined_like ${LAPACK_LIB_FLAGS})
#target_link_libraries(test_combined_like ${CLASSLIB})
#target_link_libraries(test_combined_like ${PLANCKLIB})

# compile learn_matrix
set(CMAKE_INCLUDE_CURRENT_DIR ON)
add_executable(learn_matrix learn_matrix.cpp)
target_link_libraries(learn_matrix ${COSMOPPLIB})
target_link_libraries(learn_matrix ${MPI_CXX_LIBRARIES})
target_link_libraries(learn_matrix ${LAPACK_LIB_FLAGS})

