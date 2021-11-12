cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

if (UNIX)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -msse4.1 -fPIC")
endif ()
add_library(succinct_bv STATIC bit_vector.cc naive_bit_vector.cc)
target_include_directories(succinct_bv PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
