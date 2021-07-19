#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/datastructures/SimilarityCache.hpp"
#include "../../../src/Distances.hpp"
#include "../TestUtil.hpp"
class MSITest : public testing::Test {
protected:
  std::string current_data_dir = "../test/res/data/";
  MSIMetric metric_msi;
};

/**
 * Test that MSI behaves as expected in case of identity
 */
TEST_F(MSITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  std::vector<PllSplit> vec {split};
  UniquePllMap data(vec);

  EXPECT_DOUBLE_EQ(metric_msi.evaluate(0, 0, data), phylomath::h(3,3));
  free(split());
}


/**
 * Ensure correct computation of maximum values
 */
TEST_F(MSITest, maximumtest) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  UniquePllMap map({tree1, tree2});
  const PllSplitList& splits1 = map.vectors()[0];
  const PllSplitList& splits2 = map.vectors()[1];
  tree2.alignNodeIndices(tree1);
  double result = metric_msi.maximum(splits1, splits2);
  double expected_info_content = 2 * (2 * std::log2(7) + std::log2(35.0 / 3));
}

/**
 * Test MSI on a simple constructed example
 */
TEST_F(MSITest, test_msi) {
  PllSplit::setTipCount(12);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);

  std::vector<PllSplit> vec {split_a, split_b};
  UniquePllMap data(vec);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(0, 1, data), -std::log2(1.0d/21));
  free(split_a());
  free(split_b());
}

/**
 * Test MSI on an example proposed by Luise :)
 */
TEST_F(MSITest, test_luise_graph) {
  PllSplit::setTipCount(8);
  PllSplit split_1 = TestUtil::createSplit({0, 1});
  PllSplit split_2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7});
  PllSplit split_3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7});
  PllSplit split_4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5});
  std::vector<PllSplit> vec {split_1, split_2, split_3, split_4};
  UniquePllMap data(vec);
  double result = std::log2(3);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(0, 1, data), result);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(0, 2, data), result);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(0, 3, data), result);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(1, 2, data), result);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(1, 3, data), result);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(2, 3, data), result);

  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}

/**
 * Test MSI on a minimal example (see slideshow)
 */
TEST_F(MSITest, distance_from_slideshow_msi) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[1];
  UniquePllMap map = UniquePllMap({tree1, tree2});
  SimilarityCacheLinear cache(map, metric_msi);
  const PllSplitList& s1 = map.vectors()[0];
  const PllSplitList& s2 = map.vectors()[1];
  //tree2.alignNodeIndices(tree1);

  size_t split_count = s1.getSplits().size();
  std::vector<std::vector<double>> result = std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
  Distances::similaritiesForSplits(s1, s2, cache, &result);
  double alpha = std::log2(7);
  double beta = std::log2(5);
  double gamma = std::log2(3);
  double delta = 0;
  EXPECT_DOUBLE_EQ(result[0][0], gamma);
  EXPECT_DOUBLE_EQ(result[0][1], alpha);
  EXPECT_DOUBLE_EQ(result[0][2], gamma);

  EXPECT_DOUBLE_EQ(result[1][0], delta);
  EXPECT_DOUBLE_EQ(result[1][1], beta);
  EXPECT_DOUBLE_EQ(result[1][2], beta);

  EXPECT_DOUBLE_EQ(result[2][0], gamma);
  EXPECT_DOUBLE_EQ(result[2][1], gamma);
  EXPECT_DOUBLE_EQ(result[2][2], alpha);
}
