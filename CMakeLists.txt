cmake_minimum_required(VERSION 3.17) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)

project("PowerGrid" CXX)

option(OPENACC_GPU "Enable OpenACC support on NVIDIA GPU" ON)
option(OPENACC_MP "Enable OpenACC support on CPU (PGI Only) " OFF)
option(MPISupport "Enable MPI Support" ON)

option(PythonSupport "Enable Python Support" ON)

if(OPENACC_GPU AND OPENACC_MP)
	message(FATAL_ERROR "You cannot build with both OPENACC_GPU AND OPENACC_MP at the same time")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find Boost
set(Boost_USE_STATIC_LIBS OFF)
# Boost version >= 1.70.0 have improved support for cmake built into the boost install
# Unfortunately, this support breaks when used with PGI compilers
# As a work around, ignore the generated cmake modules
#set(Boost_NO_BOOST_CMAKE ON)

if(MPISupport)
    find_package(MPI)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization mpi)
else(MPISupport)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization)
endif(MPISupport)

if(Boost_FOUND)
	set(LIBS ${LIBS} ${Boost_LIBRARIES})
	include_directories(${Boost_INCLUDE_DIRS})
endif()

#set(CUDAVER 11.5)
#set(PGI_YEAR $ENV{HPCSDK_RELEASE})
if(OPENACC_GPU) 
	find_package(CUDAToolkit)
	set(LIBS ${LIBS} CUDA::cufft)
	set(LIBS ${LIBS} CUDA::cudart)
	set(LIBS ${LIBS} CUDA::nvToolsExt)
	#find_path(cuFFT_INCLUDE_DIR cufft.h HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/math_libs/${CUDAVER}/include)
	#include_directories(${cuFFT_INCLUDE_DIR})
	#find_library(cuFFT_LIBRARIES cufft HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/math_libs/${CUDAVER}/lib64)
	#set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	#find_library(cudaRT_LIBRARIES cudart HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/cuda/${CUDAVER}/targets/x86_64-linux/lib /usr/local/cuda/lib64)
	#set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
	#find_library(nvToolsExt_LIBRARIES nvToolsExt HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/cuda/${CUDAVER}/targets/x86_64-linux/lib /usr/local/cuda/lib64)
	#set(LIBS ${LIBS} ${nvToolsExt_LIBRARIES})
	find_path(FFTW_INCLUDE_DIR fftw3.h HINTS $ENV{FFTW_INC})
	find_library(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
else()
	find_package(FFTW REQUIRED FLOAT_LIB DOUBLE_LIB)
	set(LIBS ${LIBS} ${FFTW_DOUBLE_LIB} ${FFTW_FLOAT_LIB})
	include_directories(${FFTW_INCLUDE_DIRS})
endif()

# Find Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})

# Add HDF5 libraries
find_package(HDF5 COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})

include_directories(Support/ArmaExtensions)

find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)
set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})

if(MPISupport)
	include_directories(MPI)
endif(MPISupport)

message(STATUS "Libraries to be linked in: ${LIBS}")

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)

	else()
	add_definitions(-DBOOST_NO_ALIGNMENT)
		if(OPENACC_GPU)
			add_definitions(-DOPENACC_GPU)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -cuda -Minfo=accel -fast -acc -gpu=cuda12.6,cc60,cc70,cc75,cc86,cc89 ")
		elseif(OPENACC_MP)
			add_definitions(-DOPENACC_MP)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mnollvm -O1 -Minfo=accel -acc  -gpu=multicore -pgf90libs")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -Minfo=opt -pgf90libs")
		endif()

    endif()
endif()

add_definitions(-DARMA_NO_DEBUG)
add_definitions(-DARMA_DONT_USE_HDF5)
#add_definitions(-DARMA_DONT_USE_OPENMP)
if(OPENACC_GPU)
	add_definitions(-DUSE_NVTX)
endif()

# Include cnpy
add_subdirectory(cnpy)
include_directories(${PROJECT_SOURCE_DIR}/cnpy)

add_subdirectory(PowerGrid)

add_subdirectory(apps)
