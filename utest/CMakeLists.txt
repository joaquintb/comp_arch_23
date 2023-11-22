add_executable(utest
	unit_test_1.cpp
)

add_executable(utest
		progargs_test.cpp
)

target_link_libraries(utest
	PRIVATE
	sim
	GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(utest)
