cmake_minimum_required(VERSION 2.6)
project(Terran)

include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB CoreFiles "source/*.cpp")

file(GLOB HeaderFiles "include/*.h")

set(STATIC_LIBRARY OFF CACHE BOOL "Build static as opposed to shared library")

if(STATIC_LIBRARY)
    add_definitions(-DTERRAN_BUILDING_STATIC_LIBRARY)
    add_library(Terran STATIC ${CoreFiles} ${HeaderFiles})
else(STATIC_LIBRARY)
    add_definitions(-DTERRAN_BUILDING_SHARED_LIBRARY)
    add_library(Terran SHARED ${CoreFiles} ${HeaderFiles})
endif(STATIC_LIBRARY)

add_subdirectory(tests)
