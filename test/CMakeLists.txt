# see https://cmake.org/cmake/help/latest/module/FindBoost.html
set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework iostreams)

include_directories(${Boost_INCLUDE_DIR})
include_directories("${PROJECT_SOURCE_DIR}/include")

set(HTSLIB_INCLUDE_DIR "/usr/include/htslib/")
set(HTSLIB_LIBRARY "-lhts")
include_directories(${HTSLIB_INCLUDE_DIR})

find_package(ZLIB)

# Note: a target should be already defined using 'add_executable' or 'add_library'
add_executable(Boost_Tests_run src/testMotifCount.cpp )
target_link_libraries(Boost_Tests_run ${Boost_LIBRARIES} ZLIB::ZLIB ${HTSLIB_LIBRARY})
