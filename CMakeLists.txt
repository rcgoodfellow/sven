cmake_minimum_required(VERSION 3.0)

SET(CMAKE_CXX_COMPILER clang++)
SET(CMAKE_C_COMPILER clang)

project(sven)

SET(SHARED_FLAGS "-DDEBUG -Wall -Wextra -O2 -g")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SHARED_FLAGS} -std=c++11 -fpic")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SHARED_FLAGS} -std=c11")


include_directories(
  include
  gtest/include
  /opt/intel/composer_xe_2015.0.064/tbb/include
)

link_directories(
  /opt/intel/composer_xe_2015.0.064/tbb/lib/intel64/gcc4.4
)

add_subdirectory(gtest)

add_library(sven_math SHARED src/CoreMath.cxx)
add_library(sven_runtime SHARED src/Runtime.cxx)

add_executable(OperatorTests test/OperatorTests.cxx)
target_link_libraries(OperatorTests sven_math sven_runtime tbb gtest gtest_main 
  pthread)
