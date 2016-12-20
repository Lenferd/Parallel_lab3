cmake_minimum_required(VERSION 3.6)
project(Optimization3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(MPI REQUIRED)
include_directories( ${MPI_INCLUDE_PATH} )

set(SOURCE_FILES src/main.cpp)
add_executable(Optimization3 ${SOURCE_FILES})
target_link_libraries( Optimization3 ${MPI_LIBRARIES} )
