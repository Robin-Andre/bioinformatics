#include "gtest/gtest.h"
#include "../../../src/MaximumMatcher.hpp"

#include <random>
#include <iomanip>


class MaximumMatcherTest : public testing::Test {

protected:

};


TEST_F(MaximumMatcherTest, test_assignment) {
    std::vector<std::vector<double>> weights = {{0.2d, 0.1d, 0.1d, -0.1d}, {0.0d, 0.2d, -0.1d, 0.1d}, {0.1d, 0.1d, 0.2d, 0.0d}, {0.0d, -0.1d, 0.1d, 0.2d}};
    EXPECT_DOUBLE_EQ(MaximumMatcher::match(weights), 0.8d);
}


TEST_F(MaximumMatcherTest, test_maximum) {
  size_t n = 10;
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distr(-10, 10);
  std::vector<std::vector<double>> weights = std::vector<std::vector<double>> (n, std::vector<double>(n));
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j < n; ++j){
      weights[i][j] = distr(eng);
    }
  }

  size_t assignment[weights.size()];
  for(size_t i = 0; i < weights.size(); ++i){
    assignment[i] = i;
  }
  double maximum = MaximumMatcher::match(weights);
  size_t len = sizeof(assignment) / sizeof(assignment[0]);
  do {
   double current = 0;
   for (size_t i = 0; i < n; ++i){
     current += weights[i][assignment[i]];
   }
   ASSERT_TRUE(current <= maximum);
 } while (std::next_permutation(assignment, assignment + len));

}
