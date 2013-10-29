#include_directories(${Terran_SOURCE_DIR}/include)
link_directories(${Terran_BINARY_DIR})


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Terran_BINARY_DIR})

file(GLOB test_progs "test*.cpp")
foreach(test_prog ${test_progs})
	get_filename_component(test_root ${test_prog} NAME_WE)
	add_executable(${test_root} ${test_prog} util.cpp)
	target_link_libraries(${test_root} Terran)
	add_test(${test_root} ${Terran_BINARY_DIR}/${test_root})
endforeach(test_prog ${test_progs})