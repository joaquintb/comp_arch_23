add_executable(ftest
	functional_test.cpp
)

target_link_libraries(ftest 
	PRIVATE
	sim 
	GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(ftest)
