cmake_minimum_required(VERSION 3.5)
project(MINIMIZERS)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

option(MINIMIZERS_USE_SANITIZERS "Build Debug and RelWithDebInfo with ASAN" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
message(STATUS "Compiling for processor: " ${CMAKE_HOST_SYSTEM_PROCESSOR})
if(UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  message(STATUS "Compiling with flags: -march=native -mbmi2 -msse4.2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

if(UNIX AND CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  message(STATUS "Compiling with flags: ${CMAKE_CXX_FLAGS}")
endif()

add_compile_options(-Wall -Wextra -Wno-missing-braces -Wno-unused-function)

include_directories(.)
include_directories(include)
add_library(minimizers STATIC src/algorithms.cpp)
target_include_directories(minimizers PUBLIC include)

add_executable(density tool/density.cpp)
target_compile_options(density PRIVATE -Wall -Wextra -Wno-missing-braces
                                       -Wno-unused-function)

if(MINIMIZERS_USE_SANITIZERS)
  target_compile_options(
    minimizers PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address
                       -fsanitize=undefined -fno-omit-frame-pointer>)
  target_link_options(
    minimizers PRIVATE $<$<CONFIG:Debug,RelWithDebInfo>:-fsanitize=address
    -fsanitize=undefined>)
  message(
    STATUS
      "Using sanitizers. Compiling with flags: -fsanitize=address -fno-omit-frame-pointer"
  )
endif()

target_link_libraries(density z minimizers)
add_executable(generate_random_fasta tool/generate_random_fasta.cpp)
