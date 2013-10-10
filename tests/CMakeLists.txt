#include_directories(${Terran_SOURCE_DIR}/include)
link_directories(${Terran_BINARY_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Terran_BINARY_DIR})

add_executable(testEMGaussian testEMGaussian.cpp util.cpp)
target_link_libraries(testEMGaussian Terran)

add_executable(testEMPeriodicGaussian testEMPeriodicGaussian.cpp util.cpp)
target_link_libraries(testEMPeriodicGaussian Terran)

add_executable(testMethodsPeriodicGaussian testMethodsPeriodicGaussian.cpp util.cpp)
target_link_libraries(testMethodsPeriodicGaussian Terran)

add_executable(testMethodsGaussian testMethodsGaussian.cpp util.cpp)
target_link_libraries(testMethodsGaussian Terran)

add_executable(testMathFunctions testMathFunctions.cpp util.cpp)
target_link_libraries(testMathFunctions Terran)

add_executable(testCluster testCluster.cpp util.cpp)
target_link_libraries(testCluster Terran)

add_executable(testClusterTree testClusterTree.cpp util.cpp)
target_link_libraries(testClusterTree Terran)

add_executable(testRealData testRealData.cpp util.cpp)
target_link_libraries(testRealData Terran)