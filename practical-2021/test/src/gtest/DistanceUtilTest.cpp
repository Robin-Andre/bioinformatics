#include "gtest/gtest.h"
#include "../../../src/DistanceUtil.hpp"
#include "../../../src/metrics/SPI.hpp"
#include "../../../src/metrics/MSI.hpp"
#include "../../../src/metrics/MCI.hpp"
#include "../../../src/io/TreeReader.hpp"
class DistanceUtilTest : public testing::Test {
  protected: //std::string current_test_dir = "../test/res/data/heads/BS/";
  std::string current_data_dir = "../test/res/data/";
  virtual void SetUp() {

  }
  virtual void TearDown() {

  }
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
  SPI metric_spi;
  std::vector<std::vector<double>> result = DistanceUtil::similaritiesForSplits(splits1, splits2, metric_spi);
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
  MSI metric_msi;
  std::vector<std::vector<double>> result = DistanceUtil::similaritiesForSplits(splits1, splits2, metric_msi);
  double alpha = std::log(7);
  double beta = std::log(5);
  double gamma = std::log(3);
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
  MSI metric_msi;
  double result = DistanceUtil::maximumValue(splits1, splits2, metric_msi);
  std::cout << "maximumMSI: " << result << "\n";

  SPI metric_spi;
  result = DistanceUtil::maximumValue(splits1, splits2, metric_spi);
  std::cout << "maximumSPI: " << result << "\n";

  MCI metric_mci;
  result = DistanceUtil::maximumValue(splits1, splits2, metric_mci);
  std::cout << "maximumMCI: " << result << "\n";

  
}