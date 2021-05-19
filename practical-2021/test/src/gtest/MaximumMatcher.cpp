#include "gtest/gtest.h"
#include "../../../src/MaximumMatcher.hpp"


class MaximumMatcherTest : public testing::Test {

protected:

};


TEST_F(MaximumMatcherTest, test_assignment) {
    std::vector<std::vector<double>> weights = {{2, 1, 1, -1}, {0, 2, -1, 1}, {1, 1, 2, 0}, {0, -1, 1, 2}};
    EXPECT_EQ(MaximumMatcher::match(weights), 8);
}
