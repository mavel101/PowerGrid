cmake_minimum_required(VERSION 3.19) #CMake version 3.19 required to properly support FindCUDAToolkit

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

project("PowerGrid" C CXX)

option(OPENACC_GPU "Enable OpenACC support on NVIDIA GPU" OFF)
option(OPENACC_MP "Enable OpenACC support on CPU (PGI Only) " ON)
option(MPISupport "Enable MPI Support" ON)

if(OPENACC_GPU AND OPENACC_MP)
	message(FATAL_ERROR "You cannot build with both OPENACC_GPU AND OPENACC_MP at the same time")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find Boost
set(Boost_USE_STATIC_LIBS OFF)

if(MPISupport)
    find_package(MPI)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization python38 mpi) # add Python for wrapping
else(MPISupport)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization python38)
endif(MPISupport)

if(Boost_FOUND)
	set(LIBS ${LIBS} ${Boost_LIBRARIES})
	include_directories(${Boost_INCLUDE_DIRS})
endif()

# Find python for wrapping
find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
	set(LIBS ${LIBS} ${PYTHON_LIBRARIES})
	include_directories(${Python_INCLUDE_DIRS})
endif()

set(CUDAVER 10.0)
set(PGI_YEAR 2019)

# Manage GPU libraries and FFTW
if(OPENACC_GPU)
	#find_package(CUDAToolkit)

	find_path(cuFFT_INCLUDE_DIR cufft.h HINTS $ENV{NVCOMPILERS}/$ENV{NVARCH}/$ENV{NVSDKVER}/math_libs/include)
	include_directories(${cuFFT_INCLUDE_DIR})
	find_library(cuFFT_LIBRARIES cufft HINTS $ENV{NVCOMPILERS}/$ENV{NVARCH}/$ENV{NVSDKVER}/math_libs/lib64)
	set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	find_library(cudaRT_LIBRARIES cudart HINTS $ENV{NVCOMPILERS}/$ENV{NVARCH}/$ENV{NVSDKVER}/cuda/lib64 /usr/local/cuda/lib64)
	set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
	find_library(nvToolsExt_LIBRARIES nvToolsExt HINTS $ENV{NVCOMPILERS}/$ENV{NVARCH}/$ENV{NVSDKVER}/cuda/lib64 /usr/local/cuda/lib64)
	set(LIBS ${LIBS} ${nvToolsExt_LIBRARIES})
	find_package(FFTW REQUIRED FLOAT_LIB DOUBLE_LIB)
	set(LIBS ${LIBS} ${FFTW_DOUBLE_LIB} ${FFTW_FLOAT_LIB})
	include_directories(${FFTW_INCLUDE_DIRS})
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

#Define files to be compiled
set(PG_FILES
	${PROJECT_SOURCE_DIR}/PowerGrid/Robject.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/QuadPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TVPenalty.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftCPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/fftGPU.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpu.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/ftCpuWithGrads.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gdft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/GdftR2.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gfft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/TimeSegmentation.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/griddingSupport.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/gridding.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/Gnufft.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/SENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSENSE.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/pcSenseTimeSeg.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/reconSolve.cpp
	${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.cpp
	#${PROJECT_SOURCE_DIR}/PowerGrid/acqTracking.h
	${PROJECT_SOURCE_DIR}/PowerGrid/reconSupport/directRecon.cpp
)

add_library(PowerGrid SHARED ${PG_FILES})
target_link_libraries(PowerGrid ${LIBS})

if(MPISupport)
	include_directories(MPI)
endif(MPISupport)

find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/ ~/CIVIC/lib)
set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})

message(STATUS "Libraries to be linked in: ${LIBS}")

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)

	else()
	add_definitions(-DBOOST_NO_ALIGNMENT)
		if(OPENACC_GPU)
			add_definitions(-DOPENACC_GPU)
			#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mp -Mnollvm -Mcuda -Minfo=accel -fast -acc -ta=tesla:cc35,cc60,cc75,flushz -pgf90libs")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -acc -pgf90libs")
		elseif(OPENACC_MP)
			add_definitions(-DOPENACC_MP)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mnollvm -O1 -Minfo=accel -acc  -ta=multicore -pgf90libs")
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

set(POWERGRID_ISMRMRD
	./PowerGrid/PowerGridIsmrmrd.cpp
)
add_executable(PowerGridIsmrmrd ${POWERGRID_ISMRMRD})
target_link_libraries(PowerGridIsmrmrd ${LIBS} PowerGrid)

set(POWERGRID_PCSENSE
	./PowerGrid/PowerGridPcSense.cpp
)
add_executable(PowerGridPcSense ${POWERGRID_PCSENSE})
target_link_libraries(PowerGridPcSense ${LIBS} PowerGrid)

set(POWERGRID_PCSENSE_TIMESEG
	./PowerGrid/PowerGridPcSenseTimeSeg.cpp
)
add_executable(PowerGridPcSenseTimeSeg ${POWERGRID_PCSENSE_TIMESEG})
target_link_libraries(PowerGridPcSenseTimeSeg ${LIBS} PowerGrid)

# Python bindings
find_package(pybind11)
set(POWERGRID_PY
	./PowerGrid/PowerGridPy.cpp
)
set(CMAKE_SHARED_MODULE_PREFIX "")
add_library(PowerGridPy MODULE ${POWERGRID_PY})
target_link_libraries(PowerGridPy ${LIBS} PowerGrid)
target_include_directories(PowerGrid PRIVATE ${PYTHON_INCLUDE_DIRS})

add_subdirectory(PowerGrid/lowRank)
add_subdirectory(PowerGrid/reconSupport)

if (MPISupport)
    add_subdirectory(PowerGrid/MPI)
endif (MPISupport)

# Manage installing targets to sane locations in /opt
install(TARGETS PowerGridPcSense DESTINATION bin)
install(TARGETS PowerGridIsmrmrd DESTINATION bin)
install(TARGETS PowerGrid DESTINATION lib)
