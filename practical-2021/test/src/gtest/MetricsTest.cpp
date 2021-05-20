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

  EXPECT_EQ(DistanceUtil::doublefactorial(41), 13113070457687988603440625*1.0d);
  EXPECT_EQ(DistanceUtil::doublefactorial(43), 563862029680583509947946875*1.0d);

  //EXPECT_DOUBLE_EQ(1.0d/43, (13113070457687988603440625*1.0d)/563862029680583509947946875);
}


TEST_F(MetricsTest, test_truncated) {
  EXPECT_EQ(DistanceUtil::truncatedDoublefactorial(41, 43), 43);
  EXPECT_EQ(DistanceUtil::truncatedDoublefactorial(5, 11), 7 * 9 * 11);
  EXPECT_EQ(DistanceUtil::factorialQuotient(3, 5, 9), 3.0d /(7*9));
}

TEST_F(MetricsTest, test_phylogenetic_probability) {
  PllSplit::setTipCount(8);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(2, 3), 1.0d/5);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(2, 4), 1.0d/7);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(3, 4), 1.0d/21);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(2, 2), 1.0d/3);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(3, 3), 3.0d/35);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(4, 4), 5.0d/231);
  PllSplit::setTipCount(24);
  EXPECT_DOUBLE_EQ(DistanceUtil::phylogeneticProbability(2, 22), 1.0d/43);
}

TEST_F(MetricsTest, test_msi) {
  PllSplit::setTipCount(12);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(DistanceUtil::MSI(split_a, split_b), -std::log(1.0d/21));
  free(split_a());
  free(split_b());
}


TEST_F(MetricsTest, test_shared_phylogenetic_probability) {
  PllSplit::setTipCount(6);
  EXPECT_DOUBLE_EQ(DistanceUtil::sharedPhylogeneticProbability(2, 4, 3, 3), 1.0d/35);
}

TEST_F(MetricsTest, test_spi) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0, 1};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split_a, split_b), -std::log(1.0d/7) - std::log(3.0d/35) + std::log(1.0d/35));
  free(split_a());
  free(split_b());
}

TEST_F(MetricsTest, test_entropy) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 1, 2};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_DOUBLE_EQ(DistanceUtil::entropy(split), -std::log(1.0d/2));
  free(split());
}


TEST_F(MetricsTest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(4), 1.0d/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_a, 1), 5.0d/12);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_a, 0), 7.0d/12);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_b, 1), 1.0d/2);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_b, 0), 1.0d/2);

  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_a, 1, split_b, 1), 1.0d/4);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_a, 1, split_b, 0), 1.0d/6);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_a, 0, split_b, 1), 1.0d/4);
  EXPECT_DOUBLE_EQ(DistanceUtil::clusteringProbability(split_a, 0, split_b, 0), 1.0d/3);

  double solution = ((1.0d / 4) * std::log(6.0d / 5)) + ((1.0d / 6) * std::log(4.0d / 5)) +
    ((1.0d / 4) * std::log(6.0d / 7)) + ((1.0d / 3) * std::log(8.0d / 7));
  EXPECT_NEAR(DistanceUtil::MCI(split_a, split_b), solution, 0.00000000001);
  free(split_a());
  free(split_b());
}

TEST_F(MetricsTest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_DOUBLE_EQ(DistanceUtil::MSI(split, split), DistanceUtil::h(3,3));
  EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split, split), DistanceUtil::h(3,3));
  double p_a = DistanceUtil::clusteringProbability(split, 1);
  double p_b = DistanceUtil::clusteringProbability(split, 0);
  //EXPECT_DOUBLE_EQ(DistanceUtil::MCI(split, split), p_a * std::log(6.0d/3) + p_b * std::log(6.0d/3));
  free(split());
}

TEST_F(MetricsTest, test_trivial) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2, 3, 4, 5};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<size_t> part1_c = {0, 3, 4};
  PllSplit split_c = TestUtil::createSplit(part1_c);

  EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split_a, split_c), 0);
  EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split_b, split_c), 0);
  EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split_c, split_a), 0);
  EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split_c, split_b), 0);

  EXPECT_DOUBLE_EQ(DistanceUtil::MSI(split_a, split_a), 0);
  EXPECT_DOUBLE_EQ(DistanceUtil::MSI(split_b, split_b), 0);


  free(split_a());
  free(split_b());
  free(split_c());
}


TEST_F(MetricsTest, test_special) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0, 1, 4, 5};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2, 3};
  PllSplit split_b = TestUtil::createSplit(part1_b);

  //EXPECT_DOUBLE_EQ(DistanceUtil::SPI(split_a, split_b), 0);

  free(split_a());
  free(split_b());
}
