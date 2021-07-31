# Benchmark
This folder collects benchmark runs from the test files. To enable benchmarking set 

`static constexpr bool PRINT_EXECUTION_TIME = true;` in `test/src/gtest/DistancesTest.cpp` `line 18`

If activated a file with the current `GIT_COMMIT_HASH` will be appended with the program run time to evaluate the instances of 10 trees of the dataset. This can be used to retroactively benchmark previous commits with a controlled environment. The folder is in the .gitignore file to prevent leakage of scripts and old datasets.
