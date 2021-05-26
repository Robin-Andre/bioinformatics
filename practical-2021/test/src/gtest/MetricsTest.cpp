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
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, 1), 5.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, 0), 7.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b, 1), 1.0d/2);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b, 0), 1.0d/2);

  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, 1, split_b, 1), 1.0d/4);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, 1, split_b, 0), 1.0d/6);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, 0, split_b, 1), 1.0d/4);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, 0, split_b, 0), 1.0d/3);

  double solution = ((1.0d / 4) * std::log(6.0d / 5)) + ((1.0d / 6) * std::log(4.0d / 5)) +
    ((1.0d / 4) * std::log(6.0d / 7)) + ((1.0d / 3) * std::log(8.0d / 7));
  MCI metric_mci;  
  EXPECT_NEAR(metric_mci.evaluate(split_a, split_b), solution, 0.00000000001);
  free(split_a());
  free(split_b());
}

TEST_F(MetricsTest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  MSI metric_msi;
  SPI metric_spi;
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split, split), phylomath::h(3,3));
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split, split), phylomath::h(3,3));
  double p_a = phylomath::clusteringProbability(split, 1);
  double p_b = phylomath::clusteringProbability(split, 0);
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
  MSI metric_msi;
  SPI metric_spi;
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_a, split_c), 0);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_b, split_c), 0);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_c, split_a), 0);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_c, split_b), 0);

  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split_a, split_a), 0);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split_b, split_b), 0);


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
  SPI metric_spi;
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_a, split_b), 0);

  free(split_a());
  free(split_b());
}

TEST_F(MetricsTest, test_msi) {
  PllSplit::setTipCount(12);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  MSI metric_msi;
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split_a, split_b), -std::log(1.0d/21));
  free(split_a());
  free(split_b());
}
TEST_F(MetricsTest, test_spi) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0, 1};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  SPI metric_spi;
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_a, split_b), -std::log(1.0d/7) - std::log(3.0d/35) + std::log(1.0d/35));
  free(split_a());
  free(split_b());
}
