cmake_minimum_required(VERSION 2.6)
project(Terran)

include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB CoreFiles "source/*.cpp")

set(STATIC_LIBRARY OFF CACHE BOOL "Build static as opposed to shared library")

if(STATIC_LIBRARY)
    add_definitions(-DTERRAN_BUILDING_STATIC_LIBRARY)
else(STATIC_LIBRARY)
    add_definitions(-DTERRAN_BUILDING_SHARED_LIBRARY)
endif(STATIC_LIBRARY)
add_library(Terran SHARED ${CoreFiles})
add_subdirectory(tests)
