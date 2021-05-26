#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/DistanceUtil.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
#include "../../../src/metrics/MCI.hpp"
#include "../../../src/metrics/SPI.hpp"
#include "../../../src/metrics/MSI.hpp"


class MetricsTest : public testing::Test {};



TEST_F(MetricsTest, test_double_factorial) {
  EXPECT_EQ(phylomath::doublefactorial(0), 1);
  EXPECT_EQ(phylomath::doublefactorial(1), 1);
  EXPECT_EQ(phylomath::doublefactorial(2), 2);
  EXPECT_EQ(phylomath::doublefactorial(3), 3);
  EXPECT_EQ(phylomath::doublefactorial(4), 8);
  EXPECT_EQ(phylomath::doublefactorial(5), 15);
  EXPECT_EQ(phylomath::doublefactorial(6), 48);

  EXPECT_EQ(phylomath::doublefactorial(41), 13113070457687988603440625*1.0d);
  EXPECT_EQ(phylomath::doublefactorial(43), 563862029680583509947946875*1.0d);

  //EXPECT_DOUBLE_EQ(1.0d/43, (13113070457687988603440625*1.0d)/563862029680583509947946875);
}


TEST_F(MetricsTest, test_truncated) {
  EXPECT_EQ(phylomath::truncatedDoublefactorial(41, 43), 43);
  EXPECT_EQ(phylomath::truncatedDoublefactorial(5, 11), 7 * 9 * 11);
  EXPECT_EQ(phylomath::factorialQuotient(3, 5, 9), 3.0d /(7*9));
}

TEST_F(MetricsTest, test_phylogenetic_probability) {
  PllSplit::setTipCount(8);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(2, 3), 1.0d/5);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(2, 4), 1.0d/7);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(3, 4), 1.0d/21);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(2, 2), 1.0d/3);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(3, 3), 3.0d/35);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(4, 4), 5.0d/231);
  PllSplit::setTipCount(24);
  EXPECT_DOUBLE_EQ(phylomath::phylogeneticProbability(2, 22), 1.0d/43);
}




TEST_F(MetricsTest, test_shared_phylogenetic_probability) {
  PllSplit::setTipCount(6);
  EXPECT_DOUBLE_EQ(phylomath::sharedPhylogeneticProbability(2, 4, 3, 3), 1.0d/35);
}

//The probability of a trivial split is 1, the value of h should be 0 as in -log(1) == 0 
TEST_F(MetricsTest, h_function_trivial_split) {
  PllSplit::setTipCount(4);
  PllSplit test_split = TestUtil::createSplit({0, 1, 3});
  EXPECT_EQ(phylomath::h(test_split), 0);
  //TODO free????
}

TEST_F(MetricsTest, test_entropy) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 1, 2};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_DOUBLE_EQ(phylomath::entropy(split), -std::log(1.0d/2));
  free(split()); // Why free here? Testutils is using a calloc but still... seems weird to do mem management in test cases
}

TEST_F(MetricsTest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(4), 1.0d/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_A), 5.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_B), 7.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b, Block_A), 1.0d/2);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b, Block_B), 1.0d/2);

  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_A, split_b, Block_A), 1.0d/4);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_A, split_b, Block_B), 1.0d/6);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_B, split_b, Block_A), 1.0d/4);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_B, split_b, Block_B), 1.0d/3);

  free(split_a());
  free(split_b());
}







