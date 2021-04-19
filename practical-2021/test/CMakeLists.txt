add_subdirectory(src)
add_executable(Google_Tests_run src/unittest.cpp)
target_link_libraries(Google_Tests_run gtest gtest_main)
