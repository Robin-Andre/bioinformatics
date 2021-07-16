#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../TestUtil.hpp"
class MCITest : public testing::Test {
protected:
  std::string current_data_dir = "../test/res/data/";
  MCIMetric metric_mci;
};



/**
 * Test that MCI behaves as expected in case of identity
 */
TEST_F(MCITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  std::vector<PllSplit> vec {split};
  PllPointerMap data(vec);
  double p_a = phylomath::clusteringProbability(split.partitionSizeOf(1));
  double p_b = phylomath::clusteringProbability(split.partitionSizeOf(0));
  EXPECT_DOUBLE_EQ(metric_mci.evaluate(0, 0, data), p_a * std::log2(6.0/3) + p_b * std::log2(6.0/3));
  free(split());
}


/**
 * Ensure correct computation of maximum values
 */
TEST_F(MCITest, maximumtest) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllPointerMap map({tree1, tree2});

  PllSplitList& splits1 = map.vectors()[0];
  PllSplitList& splits2 = map.vectors()[1];
  tree2.alignNodeIndices(tree1);
  double result = metric_mci.maximum(splits1, splits2);
  double expected_entropy = 2 * (2 * (std::log2(3.0) / 3 + (2.0 / 3) * std::log2(3.0 / 2)) + 1);
  EXPECT_DOUBLE_EQ(result, expected_entropy);
}


/**
 * Test MCI on a simple constructed example
 */
TEST_F(MCITest, test_mci) {
  PllSplit::setTipCount(12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(4), 1.0/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<PllSplit> vec {split_a, split_b};
  PllPointerMap data(vec);
  double solution = ((1.0 / 4) * std::log2(6.0 / 5)) + ((1.0 / 6) * std::log2(4.0 / 5)) +
    ((1.0 / 4) * std::log2(6.0 / 7)) + ((1.0 / 3) * std::log2(8.0 / 7));
  EXPECT_NEAR(metric_mci.evaluate(0, 1, data), solution, 0.00000000001);
  free(split_a());
  free(split_b());
}

/**
 * Test MCI on an example proposed by Luise :)
 */
TEST_F(MCITest, test_luise_graph) {
  PllSplit::setTipCount(8);
  PllSplit split_1 = TestUtil::createSplit({0, 1});
  PllSplit split_2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7});
  PllSplit split_3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7});
  PllSplit split_4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5});
  std::vector<PllSplit> vec {split_1, split_2, split_3, split_4};
  PllPointerMap data(vec);
  double result = 1.0 / 2 * std::log2(4.0/3) + 1.0 / 2 * std::log2(8.0 / 9); // Yeah that one was calculated by hand
  EXPECT_NEAR(metric_mci.evaluate(0, 1, data), result, 0.00000001);
  EXPECT_NEAR(metric_mci.evaluate(0, 2, data), result, 0.00000001);
  EXPECT_NEAR(metric_mci.evaluate(0, 3, data), result, 0.00000001);
  EXPECT_NEAR(metric_mci.evaluate(1, 2, data), result, 0.00000001);
  EXPECT_NEAR(metric_mci.evaluate(1, 3, data), result, 0.00000001);
  EXPECT_NEAR(metric_mci.evaluate(2, 3, data), result, 0.00000001);
  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}
