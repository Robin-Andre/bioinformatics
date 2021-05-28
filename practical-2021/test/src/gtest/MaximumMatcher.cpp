#include "gtest/gtest.h"
#include "../../../src/MaximumMatcher.hpp"
#include "../../../src/DistanceUtil.hpp"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/datastructures/PllTree.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/metrics/SPI.hpp"

#include <random>
#include <iomanip>


class MaximumMatcherTest : public testing::Test {
protected:
  std::string current_data_dir = "../test/res/data/";

  void checkAllPermutations(std::vector<std::vector<double>> weights, double maximum){
    double empirical_max = 0;
    size_t n = weights.size();
    bool found = false;
    size_t assignment[n];
    for(size_t i = 0; i < n; ++i){
      assignment[i] = i;
    }
    size_t len = sizeof(assignment) / sizeof(assignment[0]);
    do {
     double current = 0;
     for (size_t i = 0; i < n; ++i){
       current += weights[i][assignment[i]];
     }
     empirical_max = std::max(empirical_max, current);
   } while (std::next_permutation(assignment, assignment + len));
   EXPECT_DOUBLE_EQ(maximum, empirical_max);
  }

};


TEST_F(MaximumMatcherTest, test_assignment) {
    std::vector<std::vector<double>> weights = {{0.2d, 0.1d, 0.1d, -0.1d}, {0.0d, 0.2d, -0.1d, 0.1d}, {0.1d, 0.1d, 0.2d, 0.0d}, {0.0d, -0.1d, 0.1d, 0.2d}};
    EXPECT_DOUBLE_EQ(MaximumMatcher::match(weights), 0.8d);
}

//TODO move randomness out of here
TEST_F(MaximumMatcherTest, test_maximum) {
  size_t n = 10;
  std::random_device rd; 
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distr(-10, 10);
  std::vector<std::vector<double>> weights = std::vector<std::vector<double>> (n, std::vector<double>(n));
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j < n; ++j){
      weights[i][j] = distr(eng);
      //std::cout << weights[i][j] << " ";
    }
    //std::cout << "|\n";
  }
  double maximum = MaximumMatcher::match(weights);
  checkAllPermutations(weights, maximum);

}


TEST_F(MaximumMatcherTest, test_real){
  size_t n = 10;
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllSplit::setTipCount(tree.getTipCount());
  PllSplitList split_list = PllSplitList(tree);
  SPI spi_metric;
  std::vector<std::vector<double>> similarities = DistanceUtil::similaritiesForSplits(split_list, split_list, spi_metric);
  std::vector<std::vector<double>> weights = std::vector<std::vector<double>> (n, std::vector<double>(n));
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j < n; ++j){
      weights[i][j] = similarities[i][j];
    }
  }
  double maximum = MaximumMatcher::match(weights);
  checkAllPermutations(weights, maximum);
}
