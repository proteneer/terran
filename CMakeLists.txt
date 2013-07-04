cmake_minimum_required(VERSION 2.6)
project(EM)

set(source "source")

add_library(EM ${source}/EM.cpp ${source}/EMGaussian.cpp)
