#include "gtest/gtest.h"
#include "MaximumMatcher.hpp"
#include "io/TreeReader.hpp"
#include "Metric.hpp"
#include "Distances.hpp"

#include <random>
#include <iomanip>


class MaximumMatcherTest : public testing::Test {
protected:
  std::string current_data_dir = "../test/res/data/";
  MSIMetric msi;
  SPIMetric spi;
  MCIMetric mci;


  //Checks whether the calculated maximum is maximal among all possible matchings
  void checkAllPermutations(std::vector<std::vector<double>> weights, double maximum){
    double empirical_max = 0;
    size_t n = weights.size();
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
    std::vector<std::vector<double>> weights = {{0.2d, 0.1d, 0.1d, -0.1d},
                                                {0.0d, 0.2d, -0.1d, 0.1d},
                                                {0.1d, 0.1d, 0.2d, 0.0d},
                                                {0.0d, -0.1d, 0.1d, 0.2d}};
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

/**
 * Tests the matcher with real MCI weights, using the same tree twice
 */
TEST_F(MaximumMatcherTest, test_mci_weights){
  size_t n = 10;
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllSplit::setTipCount(tree.getTipCount());
  UniquePllMap map = UniquePllMap({tree});
  SimilarityCacheLinear cache(map, mci);
  const PllSplitList& split_list = map.vectors()[0];
  size_t split_count = split_list.getSplits().size();
  std::vector<std::vector<double>> similarities =
    std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
  Distances::similaritiesForSplits(split_list, split_list, cache, &similarities);
  std::vector<std::vector<double>> weights = std::vector<std::vector<double>> (n, std::vector<double>(n));
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j < n; ++j){
      weights[i][j] = similarities[i][j];
    }
  }
  double maximum = MaximumMatcher::match(weights);
  checkAllPermutations(weights, maximum);
}

/**
 * Tests the matcher with real MCI weights, using two different trees
 */
TEST_F(MaximumMatcherTest, test_mci_weights_two_trees) {
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/24")[2];
  UniquePllMap map = UniquePllMap({tree1, tree2});
  SimilarityCacheLinear cache(map, mci);
  const PllSplitList& s1 = map.vectors()[0];
  const PllSplitList& s2 = map.vectors()[1];
  PllSplit::setTipCount(tree1.getTipCount());
  size_t split_count = s1.getSplits().size();
  std::vector<std::vector<double>> similarities =
    std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
  Distances::similaritiesForSplits(s1, s2, cache, &similarities);

  //std::vector<size_t> match_results = MaximumMatcher::match_vector(similarities); We cannot trivially test matchings
  double match = MaximumMatcher::match(similarities);
  EXPECT_NEAR(match, 10.6288, 0.0001);
}
