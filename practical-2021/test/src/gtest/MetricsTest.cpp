#include "gtest/gtest.h"
#include "../../../src/PllSplits.hpp"
#include "../../../src/DistanceUtil.hpp"
#include "../TestUtil.hpp"


class MetricsTest : public testing::Test {};



TEST_F(MetricsTest, test_double_factorial) {
  EXPECT_EQ(DistanceUtil::doublefactorial(0), 1);
  EXPECT_EQ(DistanceUtil::doublefactorial(1), 1);
  EXPECT_EQ(DistanceUtil::doublefactorial(2), 2);
  EXPECT_EQ(DistanceUtil::doublefactorial(3), 3);
  EXPECT_EQ(DistanceUtil::doublefactorial(4), 8);
  EXPECT_EQ(DistanceUtil::doublefactorial(5), 15);
  EXPECT_EQ(DistanceUtil::doublefactorial(6), 48);
}

TEST_F(MetricsTest, test_phylogenetic_probability) {
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(2, 3), 1.0d/5);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(3, 4), 1.0d/21);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(2, 2), 1.0d/3);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(3, 3), 3.0d/35);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(4, 4), 5.0d/231);
}

TEST_F(MetricsTest, test_msi) {
  PllSplit::setTipCount(12);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {2, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(DistanceUtil::MSI(split_a, split_b), -std::log(1.0d/21));
  free(split_a());
  free(split_b());
}
