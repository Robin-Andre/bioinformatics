#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
class SPITest : public testing::Test {
protected:
  SPIMetric metric_spi;
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

TEST_F(SPITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  metric_spi.evaluate(test_variable, split, split);
  phylomath::h(result_variable, 3, 3);
  EXPECT_EQ(mpfr_cmp(result_variable, test_variable), 0);
  free(split());
}
/*TODO this test case only works if a phylogenetic probability of "all splits in the same partition" == 1
  essentially part1_b is causing this behaviour. Since I am unsure whether we will ever have such partitions
  I am reluctant to enable this test.

TEST_F(SPITest, test_trivial) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2, 3, 4, 5};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  std::vector<size_t> part1_c = {0, 3, 4};
  PllSplit split_c = TestUtil::createSplit(part1_c);
  SPI metric_spi;
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_a, split_c), 0);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_b, split_c), 0);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_c, split_a), 0);
  EXPECT_DOUBLE_EQ(metric_spi.evaluate(split_c, split_b), 0);

  free(split_a());
  free(split_b());
  free(split_c());
}*/

TEST_F(SPITest, test_incompatible_splits) {
  PllSplit::setTipCount(4);
  std::vector<size_t> part1_a = {0, 1};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 2};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  metric_spi.evaluate(test_variable, split_a, split_b);
  EXPECT_EQ(mpfr_sgn(test_variable), 0);
  free(split_a());
  free(split_b());

}


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
  metric_spi.evaluate(test_variable, split_a, split_b);
  double result = -2*std::log2(15.0/105)+std::log2(3.0/105);
  EXPECT_EQ(mpfr_get_d(test_variable, RND), result);
  free(split_a());
  free(split_b());
}

TEST_F(SPITest, test_spi) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1_a = {0, 1};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 1, 2};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  double result =  -std::log2(1.0d/7) - std::log2(3.0d/35) + std::log2(1.0d/35);
  metric_spi.evaluate(test_variable, split_a, split_b);
  EXPECT_EQ(mpfr_get_d(test_variable, RND), result);

  free(split_a());
  free(split_b());
}
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
  double result = -2*std::log2(1.0d/11) + std::log2(1.0d/99);
  metric_spi.evaluate(test_variable, split_1, split_2);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_spi.evaluate(test_variable, split_1, split_3);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_spi.evaluate(test_variable, split_1, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_spi.evaluate(test_variable, split_2, split_3);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_spi.evaluate(test_variable, split_2, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_spi.evaluate(test_variable, split_3, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}
