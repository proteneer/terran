include_directories(${EM_SOURCE_DIR}/source)
link_directories(${EM_BINARY_DIR}/source)

add_executable(testEMGaussian testEMGaussian.cpp util.cpp)
target_link_libraries(testEMGaussian EM)

add_executable(testEMPeriodicGaussian testEMPeriodicGaussian.cpp util.cpp)
target_link_libraries(testEMPeriodicGaussian EM)

add_executable(testPartitionPeriodicGaussian testPartitionPeriodicGaussian.cpp util.cpp)
target_link_libraries(testPartitionPeriodicGaussian EM)

add_executable(testAssignment testAssignment.cpp util.cpp)
target_link_libraries(testAssignment EM)
