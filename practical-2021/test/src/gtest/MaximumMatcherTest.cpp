#include "gtest/gtest.h"
#include "../../../src/MaximumMatcher.hpp"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/datastructures/PllTree.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/Metric.hpp"

#include <random>
#include <iomanip>


class MaximumMatcherTest : public testing::Test {
protected:
  std::string current_data_dir = "../test/res/data/";
  MSIMetric msi;
  SPIMetric spi;
  MCIMetric mci;

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
  size_t s = split_list.getSplits().size();
  std::vector<std::vector<double>>  similarities = std::vector<std::vector<double>>(s, std::vector<double>(s));
  spi.similaritiesForSplits(split_list, split_list, similarities);
  std::vector<std::vector<double>> weights = std::vector<std::vector<double>> (n, std::vector<double>(n));
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j < n; ++j){
      weights[i][j] = similarities[i][j];
    }
  }
  double maximum = MaximumMatcher::match(weights);
  checkAllPermutations(weights, maximum);
}
TEST_F(MaximumMatcherTest, test_unequal_mci) {

  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/24")[2];
  PllSplitList s1 = PllSplitList(tree1);
  PllSplitList s2 = PllSplitList(tree2);
  PllSplit::setTipCount(tree1.getTipCount());
  size_t n = s1.getSplits().size();
  std::vector<std::vector<double>>  similarities = std::vector<std::vector<double>>(n, std::vector<double>(n));
  mci.similaritiesForSplits(s1, s2, similarities);

  //std::vector<size_t> match_results = MaximumMatcher::match_vector(similarities); We cannot trivially test matchings
  double match = MaximumMatcher::match(similarities);
  EXPECT_NEAR(match, 10.6288, 0.0001);
}


//This test was for debugging the normalization values, we can first reenable it as soon as alexis speaks
//How to do the normalization.
TEST_F(MaximumMatcherTest, test_unequal_mci2) {

  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/141")[2];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/141")[5];
  PllSplit::setTipCount(tree1.getTipCount());
  PllSplitList s1 = PllSplitList(tree1);
  PllSplitList s2 = PllSplitList(tree2);
  size_t n = s1.getSplits().size();
  std::vector<std::vector<double>>  similarities = std::vector<std::vector<double>>(n, std::vector<double>(n));

  mci.similaritiesForSplits(s1, s2, similarities);
  /*for(unsigned i = 0; i < similarities.size(); ++i) {
    for(unsigned j = 0; j < similarities.size(); ++j) {
      std::cout << similarities[i][j] << " ";
    }
    std::cout << "\n";
  }*/
  std::vector<size_t> match_results = MaximumMatcher::match_vector(similarities);
  double match = MaximumMatcher::match(similarities);
  /*for(unsigned i = 0; i < match_results.size(); ++i) {
    std::cout << match_results[i] << " ";
    std::cout << "Node: " << i << " -> " << match_results[i] << "\n";
  }*/
  double maximum = mci.maximum(s1, s2);
  /*std::cout << "\nValue: " << match << "\n";
  std::cout << "Maximum: " << maximum << "\n";
  std::cout << "Fraction: " << match / maximum << "\n";
  std::cout << "Normalized(x2): " << 2*(maximum - match) << "\n";
  std::cout << "Double normalized(x2) " <<  2*(maximum - match) / (2*maximum) << "\n";*/
}
