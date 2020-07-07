cmake_minimum_required(VERSION 3.14)

project(gear LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

#temp setting for profiling
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")


include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/Modules/MotifCount")
include_directories("${PROJECT_SOURCE_DIR}/Modules/TelomereAnalysis")
include_directories("${PROJECT_SOURCE_DIR}/Modules/VariantCallAnalysis")


#EXEC SETUP
add_executable(gear src/main.cpp include/kseq.h include/utils.h include/logger.h src/utils.cpp include/csvParser.h src/Modules/MotifCount/motifCount.h src/Modules/MotifCount/motifRegion.h src/sequenceReader.cpp include/sequenceReader.h src/Modules/TelomereAnalysis/sequenceNetwork.cpp src/Modules/TelomereAnalysis/sequenceNetwork.h src/Modules/TelomereAnalysis/kmerNode.cpp src/Modules/TelomereAnalysis/kmerNode.h include/options.h src/Modules/VariantCallAnalysis/variantCallAnalysis.cpp src/Modules/VariantCallAnalysis/variantCallAnalysis.h src/Modules/MotifCount/motifCount.cpp src/Modules/MotifCount/motifRegion.cpp include/genomeRegion.h src/Modules/VariantCallAnalysis/variantRegion.cpp src/Modules/VariantCallAnalysis/variantRegion.h) #kernels.cuh ${KERNEL_O}

#BOOST SETUP
set(BOOST_ROOT "/opt/boost")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.68 COMPONENTS program_options iostreams)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()


#HSTlib SETUP
message(STATUS "Manually setting htslib up.")
if(EXISTS "/usr/include/htslib/")
    message(STATUS "Found htslib: /usr/include/htslib/")
    set(HTSLIB_INCLUDE_DIR "/usr/include/htslib/")
elseif(EXISTS "/usr/local/include/htslib/")
    message(STATUS "Found htslib: /usr/local/include/htslib/")
    set(HTSLIB_INCLUDE_DIR "/usr/local/include/htslib/")
endif()

set(HTSLIB_LIBRARY "-lhts")
include_directories(${HTSLIB_INCLUDE_DIR})

find_package(ZLIB REQUIRED)
find_package (Threads REQUIRED)

target_link_libraries(gear ${Boost_LIBRARIES} ${HTSLIB_LIBRARY} ZLIB::ZLIB ${CMAKE_THREAD_LIBS_INIT})
add_subdirectory(${PROJECT_SOURCE_DIR}/test)

