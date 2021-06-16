#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
#include <mpfr.h>
class MSITest : public testing::Test {
protected:
  MSIMetric metric_msi;
  mpfr_t test_variable;
  mpfr_t result_variable;
  virtual void SetUp() {
    mpfr_init_set_ui(test_variable, 0, RND);
    mpfr_init_set_ui(result_variable, 0, RND);
  }
  virtual void TearDown() {
    mpfr_clear(test_variable);
    mpfr_clear(result_variable);
  }
};

TEST_F(MSITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  metric_msi.evaluate(test_variable, split, split);
  phylomath::h(result_variable, 3, 3);
  EXPECT_EQ(mpfr_cmp(result_variable, test_variable), 0);
  free(split());
}
/*TODO Same as SPITest, split_b is all taxa in one partition. Is the probability of that occurence
  1 or 0. I don't know but as long as there is no definitive answer this test will stay disabled
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
}*/
TEST_F(MSITest, test_msi) {
  PllSplit::setTipCount(12);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  metric_msi.evaluate(test_variable, split_a, split_b);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), -std::log2(1.0/21));
  //EXPECT_EQ(mpfr_cmp(result_variable, test_variable), 0);
  free(split_a());
  free(split_b());
}
TEST_F(MSITest, test_luise_graph) {
  PllSplit::setTipCount(8);
  PllSplit split_1 = TestUtil::createSplit({0, 1});
  PllSplit split_2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7});
  PllSplit split_3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7});
  PllSplit split_4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5});
  MSIMetric metric_msi;
  metric_msi.evaluate(test_variable, split_1, split_2);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), std::log2(3));
  metric_msi.evaluate(test_variable, split_1, split_3);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), std::log2(3));
  metric_msi.evaluate(test_variable, split_1, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), std::log2(3));
  metric_msi.evaluate(test_variable, split_2, split_3);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), std::log2(3));
  metric_msi.evaluate(test_variable, split_2, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), std::log2(3));
  metric_msi.evaluate(test_variable, split_3, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), std::log2(3));
  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}
