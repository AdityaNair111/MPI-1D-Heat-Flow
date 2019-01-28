cmake_minimum_required(VERSION 3.9)
set(CMAKE_CXX_STANDARD 17)

project(heat1D)
include_directories(${MPI_INCLUDE_PATH})



add_executable(heat1D src/heat1D.cc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")


find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
target_link_libraries(heat1D ${MPI_C_LIBRARIES}) 
