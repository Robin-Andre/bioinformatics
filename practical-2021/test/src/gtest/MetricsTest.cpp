#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/datastructures/PllPointerMap.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/datastructures/IntersectionCache.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/Solver.hpp"
#include "../TestUtil.hpp"
#include <gmp.h>

class MetricsTest : public testing::Test {
  protected:


  std::string current_data_dir = "../test/res/data/";
  MSIMetric msi;
  SPIMetric spi;
  MCIMetric mci;

  double evaluation(double h1, double h2, double h_shared) {
      return h1 + h2 - h_shared;
  }
};

TEST_F(MetricsTest, distances_example_from_slideshow_spi) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[1];
  tree2.alignNodeIndices(tree1);
  PllPointerMap map({tree1, tree2});
  IntersectionCacheLinear cache(map, spi);
  std::vector<PllSplitList>& vec = map.vectors();
  PllSplitList& splits1 = vec[0];
  PllSplitList& splits2 = vec[1];

  size_t split_count = splits1.getSplits().size();
  std::vector<std::vector<double>> result = std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
  Solver::similaritiesForSplits(splits1, splits2, &result, cache);
  double h_standard = phylomath::h(2, 4);
  double h_i1 = phylomath::h(3, 3);
  double h_shared_beta = phylomath::h(3, 2, 6);
  double h_shared_gamma = phylomath::h(2, 2, 6);
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
TEST_F(MetricsTest, distance_from_slideshow_msi) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[1];
  PllPointerMap map = PllPointerMap({tree1, tree2});
  IntersectionCacheLinear cache(map, msi);
  PllSplitList& s1 = map.vectors()[0];
  PllSplitList& s2 = map.vectors()[1];
  //tree2.alignNodeIndices(tree1);

  size_t split_count = s1.getSplits().size();
  std::vector<std::vector<double>> result = std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
  Solver::similaritiesForSplits(s1, s2, &result, cache);
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
TEST_F(MetricsTest, maximumtest) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllPointerMap map({tree1, tree2});

  PllSplitList& splits1 = map.vectors()[0];
  PllSplitList& splits2 = map.vectors()[1];
  tree2.alignNodeIndices(tree1);
  double result = msi.maximum(splits1, splits2);
  double expected_info_content = 2 * (2 * std::log2(7) + std::log2(35.0 / 3));
  EXPECT_DOUBLE_EQ(result, expected_info_content);
  result = spi.maximum(splits1, splits2);
  EXPECT_DOUBLE_EQ(result, expected_info_content);
  result = mci.maximum(splits1, splits2);
  double expected_entropy = 2 * (2 * (std::log2(3.0) / 3 + (2.0 / 3) * std::log2(3.0 / 2)) + 1);
  EXPECT_DOUBLE_EQ(result, expected_entropy);
}


/*TEST_F(MetricsTest, rf_test) {
  PllSplit::setTipCount(64);
  std::vector<std::vector<size_t>> fst_part1s = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 1, 2, 3, 4, 5},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3}
            };
  PllSplitList fst_splitlist = TestUtil::createSplitList(fst_part1s);

  std::vector<std::vector<size_t>> snd_part1s = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList snd_splitlist = TestUtil::createSplitList(snd_part1s);

  std::vector<std::vector<size_t>> trd_part1s = {
                { 0, 2, 3, 4, 5, 6, 7},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList trd_splitlist = TestUtil::createSplitList(trd_part1s);

  RFMetric rf;
  ASSERT_EQ(rf.distanceOf(fst_splitlist, snd_splitlist, ABSOLUTE), 6);
  ASSERT_EQ(rf.distanceOf(snd_splitlist, fst_splitlist, ABSOLUTE), 6);

  ASSERT_EQ(rf.distanceOf(fst_splitlist, trd_splitlist, ABSOLUTE), 9);
  ASSERT_EQ(rf.distanceOf(trd_splitlist, fst_splitlist, ABSOLUTE), 9);

  ASSERT_EQ(rf.distanceOf(snd_splitlist, trd_splitlist, ABSOLUTE), 3);
  ASSERT_EQ(rf.distanceOf(trd_splitlist, snd_splitlist, ABSOLUTE), 3);

}*/
