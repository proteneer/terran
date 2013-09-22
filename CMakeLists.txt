cmake_minimum_required(VERSION 2.6)
project(Terran)

include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB CoreFiles "source/*.cpp")
add_library(Terran ${CoreFiles})
add_subdirectory(tests)
