# see https://cmake.org/cmake/help/latest/module/FindBoost.html
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework iostreams)

include_directories(${Boost_INCLUDE_DIR})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/Modules/MotifCount")
include_directories("${PROJECT_SOURCE_DIR}/Modules/TelomereAnalysis")
include_directories("${PROJECT_SOURCE_DIR}/Modules/VariantCallAnalysis")

#HSTlib SETUP
message(STATUS "Manually setting htslib up.")
if(EXISTS "/usr/include/htslib/")
    message(STATUS "Found htslib: /usr/include/htslib/")
    set(HTSLIB_INCLUDE_DIR "/usr/include/htslib/")
    set(HTSLIB_LIBRARY "-lhts")
elseif(EXISTS "/usr/local/include/htslib/")
    message(STATUS "Found htslib: /usr/local/include/htslib/")
    set(HTSLIB_INCLUDE_DIR "/usr/local/include/htslib/")
    set(HTSLIB_LIBRARY "-lhts")
elseif(EXISTS "/sw/include/htslib/")
    message(STATUS "Found htslib: /sw/include/htslib/")
    set(HTSLIB_INCLUDE_DIR "/sw/include/htslib/")
    set(HTSLIB_LIBRARY "-L/sw/lib/ -lhts")
endif()


include_directories(${HTSLIB_INCLUDE_DIR})

find_package(ZLIB)
find_package (Threads REQUIRED)

configure_file(${PROJECT_SOURCE_DIR}/test/data/par/input.json ${CMAKE_CURRENT_BINARY_DIR}/data/input.json COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/par/transition_matrix.csv ${CMAKE_CURRENT_BINARY_DIR}/data/transition_matrix.csv COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/input/input.fq ${CMAKE_CURRENT_BINARY_DIR}/data/input.fq COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/input/input.fq.gz ${CMAKE_CURRENT_BINARY_DIR}/data/input.fq.gz COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/input/input.csv ${CMAKE_CURRENT_BINARY_DIR}/data/input.csv COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/input/input.bam ${CMAKE_CURRENT_BINARY_DIR}/data/input.bam COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/input/input.bam.bai ${CMAKE_CURRENT_BINARY_DIR}/data/input.bam.bai COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/output/result.json ${CMAKE_CURRENT_BINARY_DIR}/data/result.json COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/test/data/output/bam_result.json ${CMAKE_CURRENT_BINARY_DIR}/data/bam_result.json COPYONLY)


# Note: a target should be already defined using 'add_executable' or 'add_library'
add_executable(Boost_Tests_run
        src/testMotifCount.cpp
        src/testUtils.cpp
        src/testSequenceReader.cpp
        ${PROJECT_SOURCE_DIR}/src/sequenceReader.cpp
        ${PROJECT_SOURCE_DIR}/src/Modules/GenomeAnalysis/sequenceNetwork.cpp
        ${PROJECT_SOURCE_DIR}/src/Modules/GenomeAnalysis/sequenceNetwork.h
        ${PROJECT_SOURCE_DIR}/src/Modules/GenomeAnalysis/kmerNode.cpp
        ${PROJECT_SOURCE_DIR}/src/Modules/GenomeAnalysis/kmerNode.h
        ${PROJECT_SOURCE_DIR}/src/Modules/MotifCount/motifCount.cpp
        ${PROJECT_SOURCE_DIR}/src/Modules/MotifCount/motifRegion.cpp
        ${PROJECT_SOURCE_DIR}/src/Modules/VariantCallAnalysis/variantRegion.cpp
        ${PROJECT_SOURCE_DIR}/src/Modules/VariantCallAnalysis/variantRegion.h
        src/testKmerNode.cpp
        src/testSequenceNetwork.cpp
        src/testGenomeRegion.cpp
        src/testVariantRegion.cpp
        src/testMotifRegion.cpp)

target_link_libraries(Boost_Tests_run ${Boost_LIBRARIES} ZLIB::ZLIB ${HTSLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
