add_executable(utest
	simulation_test.cpp
	progargs_test.cpp
	grid_test.cpp
)

target_link_libraries(utest
	PRIVATE
	sim
	GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(utest)
