#include "gtest/gtest.h"
#include "../../../src/Metric.hpp"
#include "../../../src/io/TreeReader.hpp"
class DistanceUtilTest : public testing::Test {
  protected: //std::string current_test_dir = "../test/res/data/heads/BS/";
  std::string current_data_dir = "../test/res/data/";
  MSIMetric msi;
  SPIMetric spi;
  MCIMetric mci;
  //@Softwipe, not used so not enabled
  /*virtual void SetUp() {

  }
  virtual void TearDown() {

  }*/
  double evaluation(double h1, double h2, double h_shared) {
      return h1 + h2 - h_shared;
  }
};
TEST_F(DistanceUtilTest, distances_example_from_slideshow_spi) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[1];
  PllSplitList splits1 = PllSplitList(tree1);
  PllSplitList splits2 = PllSplitList(tree2);
  tree2.alignNodeIndices(tree1);
  std::vector<std::vector<double>> result = spi.similaritiesForSplits(splits1, splits2);
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
TEST_F(DistanceUtilTest, distance_from_slideshow_msi) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[1];
  PllSplitList splits1 = PllSplitList(tree1);
  PllSplitList splits2 = PllSplitList(tree2);
  tree2.alignNodeIndices(tree1);
  std::vector<std::vector<double>> result = msi.similaritiesForSplits(splits1, splits2);
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
TEST_F(DistanceUtilTest, maximumtest) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllSplitList splits1 = PllSplitList(tree1);
  PllSplitList splits2 = PllSplitList(tree2);
  tree2.alignNodeIndices(tree1);
  double result = msi.maximumValue(splits1, splits2);
  double expected_info_content = 2 * std::log2(7) + std::log2(35.0 / 3);
  EXPECT_DOUBLE_EQ(result, expected_info_content);
  result = spi.maximumValue(splits1, splits2);
  EXPECT_DOUBLE_EQ(result, expected_info_content);
  result = mci.maximumValue(splits1, splits2);
  double expected_entropy = 2 * (std::log2(3.0) / 3 + (2.0 / 3) * std::log2(3.0 / 2)) + 1;
  EXPECT_DOUBLE_EQ(result, expected_entropy);
}
