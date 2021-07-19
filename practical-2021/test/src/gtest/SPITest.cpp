#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/datastructures/SimilarityCache.hpp"
#include "../../../src/Distances.hpp"
#include "../TestUtil.hpp"
class SPITest : public testing::Test {
protected:
  std::string current_data_dir = "../test/res/data/";
  SPIMetric metric_spi;

  double evaluation(double h1, double h2, double h_shared) {
      return h1 + h2 - h_shared;
  }
};


/**
 * Test that SPI behaves as expected in case of identity
 */
TEST_F(SPITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  std::vector<PllSplit> vec {split};
  UniquePllMap map(vec);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 0, map), phylomath::h(3,3));
  free(split());
}

/**
 * Ensure correct computation of maximum values
 */
TEST_F(SPITest, maximumtest) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  UniquePllMap map({tree1, tree2});

  PllSplitList& splits1 = map.vectors()[0];
  PllSplitList& splits2 = map.vectors()[1];
  tree2.alignNodeIndices(tree1);
  double expected_info_content = 2 * (2 * std::log2(7) + std::log2(35.0 / 3));
  double result = metric_spi.maximum(splits1, splits2);
  EXPECT_DOUBLE_EQ(result, expected_info_content);
}


/**
 * Test SPI on a simple constructed example
 */
TEST_F(SPITest, test_spi) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0, 1};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<PllSplit> vec {split_a, split_b};
  UniquePllMap map(vec);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 1, map), -std::log2(1.0d/7) - std::log2(3.0d/35) + std::log2(1.0d/35));
  free(split_a());
  free(split_b());
}

/**
 * Test SPI on an example proposed by Luise :)
 */
/*
  The following probabilities are expected:
  P(part_a) == P(part_b) = 1/11
  P(intersect) = 1/99
*/
TEST_F(SPITest, test_luise_graph) {
  PllSplit::setTipCount(8);
  PllSplit split_1 = TestUtil::createSplit({0, 1});
  PllSplit split_2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7});
  PllSplit split_3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7});
  PllSplit split_4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5});
  std::vector<PllSplit> vec {split_1, split_2, split_3, split_4};
  UniquePllMap data(vec);
  double expected_probability_single = 1.0 / 11;
  double expected_probability_intersect = 1.0 / 99;
  double result = -2 * std::log2(expected_probability_single) + std::log2(expected_probability_intersect);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 1, data), result);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 2, data), result);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 3, data), result);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(1, 2, data), result);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(1, 3, data), result);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(2, 3, data), result);
  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}

/**
 * Test SPI on a minimal example (see slideshow)
 */
TEST_F(SPITest, example_from_slideshow) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[1];
  tree2.alignNodeIndices(tree1);
  UniquePllMap map({tree1, tree2});
  SimilarityCacheLinear cache(map, metric_spi);
  std::vector<PllSplitList>& vec = map.vectors();
  PllSplitList& splits1 = vec[0];
  PllSplitList& splits2 = vec[1];

  size_t split_count = splits1.getSplits().size();
  std::vector<std::vector<double>> result = std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
  Distances::similaritiesForSplits(splits1, splits2, cache, &result);
  double h_standard = phylomath::h(2, 4);
  double h_i1 = phylomath::h(3, 3);
  double h_shared_beta = phylomath::h_shared(3, 2);
  double h_shared_gamma = phylomath::h_shared(2, 2);
  EXPECT_EQ(result[0][0], evaluation(h_standard, h_standard, h_shared_gamma));
  EXPECT_EQ(result[0][1], h_standard);
  EXPECT_EQ(result[0][2], evaluation(h_standard, h_standard, h_shared_gamma));

  EXPECT_EQ(result[1][0], 0);
  EXPECT_EQ(result[1][1], evaluation(h_i1, h_standard, h_shared_beta));
  EXPECT_EQ(result[1][2], evaluation(h_i1, h_standard, h_shared_beta));

  EXPECT_EQ(result[2][0], evaluation(h_standard, h_standard, h_shared_gamma));
  EXPECT_EQ(result[2][1], evaluation(h_standard, h_standard, h_shared_gamma));
  EXPECT_EQ(result[2][2], h_standard);
}


/**
 * Test SPI on a simple constructed example with incompatible splits
 */
TEST_F(SPITest, test_incompatible_splits) {
  PllSplit::setTipCount(4);
  std::vector<size_t> part1_a = {0, 1};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 2};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<PllSplit> vec {split_a, split_b};
  UniquePllMap map(vec);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 1 , map), 0);
  free(split_a());
  free(split_b());

}


/**
 * Test SPI for special constructed example
 */
/*Edit -> These splits are compatible: The following values are expected
  P(part_a) = 15/105
  P(part_b) =
  P(part_a & part_b) = 3 / 105
  */
TEST_F(SPITest, test_special) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0, 1, 4, 5};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2, 3};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<PllSplit> vec {split_a, split_b};
  UniquePllMap map(vec);
  double h_a = -std::log2(15.0/105);
  double h_b = -std::log2(15.0/105);
  double h_a_intersect_b = -std::log2(3.0/105);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(0, 1, map),h_a + h_b - h_a_intersect_b);
  free(split_a());
  free(split_b());
}
