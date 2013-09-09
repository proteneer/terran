include_directories(${Terran_SOURCE_DIR}/source)
link_directories(${Terran_BINARY_DIR}/source)

add_executable(testEMGaussian testEMGaussian.cpp util.cpp)
target_link_libraries(testEMGaussian EM)

add_executable(testEMPeriodicGaussian testEMPeriodicGaussian.cpp util.cpp)
target_link_libraries(testEMPeriodicGaussian EM)

add_executable(testMethodsPeriodicGaussian testMethodsPeriodicGaussian.cpp util.cpp)
target_link_libraries(testMethodsPeriodicGaussian EM)

add_executable(testMethodsGaussian testMethodsGaussian.cpp util.cpp)
target_link_libraries(testMethodsGaussian EM)

add_executable(testMathFunctions testMathFunctions.cpp util.cpp)
target_link_libraries(testMathFunctions EM)

add_executable(testCluster testCluster.cpp util.cpp)
target_link_libraries(testCluster EM)

add_executable(testClusterTree testClusterTree.cpp util.cpp)
target_link_libraries(testClusterTree EM)

add_executable(testMesh testMesh.cpp util.cpp)
target_link_libraries(testMesh EM)
