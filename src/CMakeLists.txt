cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

add_library(succinct_bv STATIC bit_vector.cc naive_bit_vector.cc)
target_include_directories(succinct_bv PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
