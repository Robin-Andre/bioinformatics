#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
class MCITest : public testing::Test {};


//This is a testcase from the old Metricstest that contained MCI, it is here now but needs a better name
TEST_F(MCITest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(4), 1.0/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  double solution = ((1.0 / 4) * std::log2(6.0 / 5)) + ((1.0 / 6) * std::log2(4.0 / 5)) +
    ((1.0 / 4) * std::log2(6.0 / 7)) + ((1.0 / 3) * std::log2(8.0 / 7));
  MCI metric_mci;  
  EXPECT_NEAR(metric_mci.evaluate(split_a, split_b), solution, 0.00000000001);
  free(split_a());
  free(split_b());
}
TEST_F(MCITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  MCI metric_mci;
  double p_a = phylomath::clusteringProbability(split, Block_A);
  double p_b = phylomath::clusteringProbability(split, Block_B);
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split, split), p_a * std::log2(6.0/3) + p_b * std::log2(6.0/3));
  free(split());
}
TEST_F(MCITest, test_luise_graph) {
  PllSplit::setTipCount(8);
  PllSplit split_1 = TestUtil::createSplit({0, 1});
  PllSplit split_2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7});
  PllSplit split_3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7});
  PllSplit split_4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5});
  MCI metric_mci;
  double result = 1.0 / 2 * std::log2(4.0/3) + 1.0 / 2 * std::log2(8.0 / 9); // Yeah that one was calculated by hand
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split_1, split_2), result);
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split_1, split_3), result);
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split_1, split_4), result);
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split_2, split_3), result);
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split_2, split_4), result);
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(split_3, split_4), result);
  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}

