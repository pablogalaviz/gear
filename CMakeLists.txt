cmake_minimum_required(VERSION 3.15)

set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

project(gear LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 14)

#temp setting for profiling
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

include_directories("${PROJECT_SOURCE_DIR}/include")

#CUDA SETUP
include(CheckLanguage)
check_language(CUDA)
find_package(CUDA  REQUIRED)
include_directories("${CUDA_INCLUDE_DIRS}")
#cuda_compile(KERNEL_O kernels.cu)

include_directories(${HDF5_CXX_INCLUDE_DIRS})


#EXEC SETUP
add_executable(gear src/main.cpp include/kseq.h include/utils.h include/logger.h src/utils.cpp include/csvParser.h include/motifCount.h include/motif.h) #kernels.cuh ${KERNEL_O}

#BOOST SETUP
set(BOOST_ROOT "/opt/boost")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.68 COMPONENTS program_options iostreams)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()


#MPI SETUP
find_package(MPI 3.1 REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

if(MPI_COMPILE_FLAGS)
    set_target_properties(gear PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(gear PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#OPENMP SETUP
find_package(OpenMP REQUIRED)

#HDF5 SETUP
find_package(HDF5 REQUIRED COMPONENTS CXX)
set(HDF5_PREFER_PARALLEL ON)
include_directories(${HDF5_CXX_INCLUDE_DIRS})

#HSTlib SETUP
message(Warning " Manually setting htslib up.")
set(HTSLIB_INCLUDE_DIR "/usr/include/htslib/")
set(HTSLIB_LIBRARY "-lhts")

include_directories(${HTSLIB_INCLUDE_DIR})


#message(STATUS "HDF5_CXX_LIBRARIES ${HDF5_CXX_LIBRARIES} ")

target_link_libraries(gear ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} OpenMP::OpenMP_CXX ${HDF5_CXX_LIBRARIES} ${HTSLIB_LIBRARY})
add_subdirectory(${PROJECT_SOURCE_DIR}/test)