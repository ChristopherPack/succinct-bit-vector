cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

if(NOT TARGET gtest)
set(GTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/lib/googletest/googletest/)
include_directories(${GTEST_ROOT}/include ${CMAKE_CURRENT_SOURCE_DIR}/../include)
link_directories(${GTEST_ROOT})
add_subdirectory(${GTEST_ROOT})
endif()

add_executable(test_naive_bit_vector
  ${CMAKE_CURRENT_SOURCE_DIR}/test_naive_bit_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/naive_bit_vector.cc)
if(UNIX)
    target_link_libraries(test_naive_bit_vector gtest gtest_main pthread)
else()
target_link_libraries(test_naive_bit_vector gtest gtest_main)
endif()

add_executable(test_bit_vector
  ${CMAKE_CURRENT_SOURCE_DIR}/test_bit_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/bit_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/naive_bit_vector.cc)
if(UNIX)
target_link_libraries(test_bit_vector gtest gtest_main pthread)
else()
target_link_libraries(test_bit_vector gtest gtest_main)
endif()