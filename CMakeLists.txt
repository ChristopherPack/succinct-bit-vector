cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(succinct_bv)

set (CMAKE_CXX_STANDARD 17)

if(UNIX)
set (CMAKE_CXX_FLAGS "-Wall -O3 -DNDEBUG -march=native -msse4.2 -mpopcnt")
else()
set (CMAKE_CXX_FLAGS_DEBUG "/MTd")
set (CMAKE_CXX_FLAGS_RELEASE "/MT")
endif()

add_subdirectory(src)
add_subdirectory(tests)
enable_testing()

add_test(NAME NaiveBitVectorTest COMMAND test_naive_bit_vector)
add_test(NAME BitVectorTest COMMAND test_bit_vector)

