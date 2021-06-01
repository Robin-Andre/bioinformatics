#include "gtest/gtest.h"
#include "../../../src/MaximumMatcher.hpp"
#include "../../../src/DistanceUtil.hpp"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/datastructures/PllTree.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/metrics/SPI.hpp"
<<<<<<< 245ada67f63b151582ee3e3808f98342d200bc47
#include "../../../src/metrics/MCI.hpp"
=======
>>>>>>> Add permutation test
#include "../../../src/metrics/MSI.hpp"

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
<<<<<<< 245ada67f63b151582ee3e3808f98342d200bc47
TEST_F(MaximumMatcherTest, test_unequal_mci) {

  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/24")[2];
  PllSplitList s1 = PllSplitList(tree1);
  PllSplitList s2 = PllSplitList(tree2);
  MCI metric_mci;
  PllSplit::setTipCount(tree1.getTipCount());
  std::vector<std::vector<double>> similarities = DistanceUtil::similaritiesForSplits(s1, s2, metric_mci);

  std::vector<size_t> match_results = MaximumMatcher::match_vector(similarities);
  double match = MaximumMatcher::match(similarities);
  for(unsigned i = 0; i < match_results.size(); ++i) {
    std::cout << match_results[i] << " ";
    //std::cout << "Node: " << i << " -> " << match_results[i] << "\n";
  }
  std::cout << "\nValue: " << match << "\n";
}
TEST_F(MaximumMatcherTest, test_unequal_msi) {

  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/24")[2];
  PllSplitList s1 = PllSplitList(tree1);
  PllSplitList s2 = PllSplitList(tree2);
  MSI metric_msi;
  PllSplit::setTipCount(tree1.getTipCount());
  std::vector<std::vector<double>> similarities = DistanceUtil::similaritiesForSplits(s1, s2, metric_msi);
  for(unsigned i = 0; i < similarities.size(); ++i) {
    for(unsigned j = 0; j < similarities.size(); ++j) {
      std::cout << similarities[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::vector<size_t> match_results = MaximumMatcher::match_vector(similarities);
  double match = MaximumMatcher::match(similarities);
  for(unsigned i = 0; i < match_results.size(); ++i) {
    std::cout << match_results[i] << " ";
    //std::cout << "Node: " << i << " -> " << match_results[i] << "\n";
  }
  double maximum = DistanceUtil::maximumValue(s1, s2, metric_msi);
  std::cout << "\nValue: " << match << "\n";
  std::cout << "Maximum: " << maximum << "\n";
  std::cout << "Fraction: " << match / maximum << "\n";
  std::cout << "Normalized(x2): " << 2*(maximum - match) << "\n";
  std::cout << "Double normalized(x2) " <<  2*(maximum - match) / (2*maximum) << "\n";
=======

TEST_F(MaximumMatcherTest, test_against_sequence){
  std::vector<size_t> reference_permutation = {0, 1, 2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 4,  5,  6,  7,  8, 9, 10, 19, 20};
  std::vector<PllTree> trees = TreeReader::readTreeFile(current_data_dir + "heads/24");
  PllTree tree = trees[0];
  PllTree tree2 = trees[2];
  PllSplit::setTipCount(tree.getTipCount());
  PllSplitList split_list = PllSplitList(tree);
  PllSplitList split_list2 = PllSplitList(tree2);
  MSI msi_metric;
  std::vector<std::vector<double>> weights = DistanceUtil::similaritiesForSplits(split_list, split_list2, msi_metric);
  std::vector<size_t> permutation = MaximumMatcher::matchingPermutation(weights);
  double maximum = MaximumMatcher::match(weights);
  double reference_weight = 0;
  for (size_t i = 0; i < weights.size(); ++i){
    reference_weight += weights[i][reference_permutation[i]];
  }
  EXPECT_EQ(reference_weight, maximum);
  EXPECT_EQ(permutation, reference_permutation);

>>>>>>> Add permutation test
}
