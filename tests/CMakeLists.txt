include_directories(${EM_SOURCE_DIR}/source)
link_directories(${EM_BINARY_DIR}/source)
add_executable(testEMGaussian testEMGaussian.cpp)
target_link_libraries(testEMGaussian EM)
