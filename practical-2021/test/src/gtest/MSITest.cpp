#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/metrics/MSI.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
class MSITest : public testing::Test {};

TEST_F(MSITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  MSI metric_msi;
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split, split), phylomath::h(3,3));
  free(split());
}

TEST_F(MSITest, test_trivial) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2, 3, 4, 5};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<size_t> part1_c = {0, 3, 4};
  PllSplit split_c = TestUtil::createSplit(part1_c);
  MSI metric_msi;
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split_a, split_a), 0);
  EXPECT_DOUBLE_EQ(metric_msi.evaluate(split_b, split_b), 0);


  free(split_a());
  free(split_b());
  free(split_c());
}
TEST_F(MSITest, test_msi) {
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