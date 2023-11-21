add_executable(utest
	unit_test_1.cpp
)

# Andrew NOTE: Can remove this when building on avagion. I only added this because I needed to include gtests locally 
include(FetchContent)
FetchContent_Declare(
  	googletest
  	URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

target_link_libraries(utest
	PRIVATE
	sim
	GTest::gtest_main
	# Microsoft.GSL::GSL
)

include(GoogleTest)
gtest_discover_tests(utest)

# 
# include(FetchContent)
# FetchContent_Declare(
#   	googletest
#   	URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
# )
# FetchContent_MakeAvailable(googletest)
# 
# enable_testing()
# 
# add_executable(
# 	unit_test_1
# 	unit_test_1.cpp
# )
# 
# target_link_libraries(
#   	unit_test_1
#  	GTest::gtest_main
# )
# include(GoogleTest)
# 