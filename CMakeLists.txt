cmake_minimum_required(VERSION 2.6)
project(Terran)

include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB CoreFiles "source/*.cpp")

file(GLOB HeaderFiles "include/*.h")

set(STATIC_LIBRARY OFF CACHE BOOL "Build static as opposed to shared library")

if(WIN32)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /openmp" )
endif(WIN32)

if(APPLE)
set(CMAKE_CXX_COMPILER "/usr/local/opt/gcc49/bin/g++-4.9")
endif(APPLE)

if(UNIX)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fopenmp" )
endif(UNIX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Valid Options: Debug, Release, RelWithDebIfno" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(STATIC_LIBRARY)
    add_definitions(-DTERRAN_BUILDING_STATIC_LIBRARY)
    add_library(Terran STATIC ${CoreFiles} ${HeaderFiles})
else(STATIC_LIBRARY)
    add_definitions(-DTERRAN_BUILDING_SHARED_LIBRARY)
    add_library(Terran SHARED ${CoreFiles} ${HeaderFiles})
endif(STATIC_LIBRARY)

enable_testing()

add_subdirectory(tests)
