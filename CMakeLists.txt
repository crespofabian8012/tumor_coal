# define the project name
project(tumor_coal)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_CC_COMPILER /usr/local/gcc-9.1/bin/gcc-9.1)
set(CMAKE_CXX_COMPILER /usr/local/gcc-9.1/bin/g++-9.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set( CMAKE_C_FLAGS "-Wall -g" )

enable_testing()

#add_subdirectory(inference)
add_subdirectory(simul)
