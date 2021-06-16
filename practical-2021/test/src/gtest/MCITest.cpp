#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
#include <mpfr.h>
class MCITest : public testing::Test {
protected:
  MCIMetric metric_mci;
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


//This is a testcase from the old Metricstest that contained MCI, it is here now but needs a better name
TEST_F(MCITest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  phylomath::clusteringProbability(test_variable, 4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  double result = ((1.0 / 4) * std::log2(6.0 / 5)) + ((1.0 / 6) * std::log2(4.0 / 5)) +
    ((1.0 / 4) * std::log2(6.0 / 7)) + ((1.0 / 3) * std::log2(8.0 / 7));
  metric_mci.evaluate(test_variable, split_a, split_b);
  EXPECT_NEAR(mpfr_get_d(test_variable, RND), result, 0.00000001);
  free(split_a());
  free(split_b());
}
TEST_F(MCITest, test_identity) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 3, 4};
  PllSplit split = TestUtil::createSplit(part1);
  mpfr_t p_a, p_b;
  mpfr_init_set_ui(p_a, 0, RND);
  mpfr_init_set_ui(p_b, 0, RND);
  phylomath::clusteringProbability(p_a, split, Block_A);
  phylomath::clusteringProbability(p_b, split, Block_B);
  mpfr_add(result_variable, p_a, p_b, RND);
  mpfr_clear(p_a);
  mpfr_clear(p_b);
  metric_mci.evaluate(test_variable, split, split);
  EXPECT_EQ(mpfr_cmp(result_variable, test_variable), 0);
  //EXPECT_DOUBLE_EQ(metric_mci.evaluate(split, split), p_a * std::log2(6.0/3) + p_b * std::log2(6.0/3));
  free(split());
}
TEST_F(MCITest, test_luise_graph) {
  PllSplit::setTipCount(8);
  PllSplit split_1 = TestUtil::createSplit({0, 1});
  PllSplit split_2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7});
  PllSplit split_3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7});
  PllSplit split_4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5});
  MCIMetric metric_mci;

  double result = 1.0 / 2 * std::log2(4.0/3) + 1.0 / 2 * std::log2(8.0 / 9); // Yeah that one was calculated by hand
  metric_mci.evaluate(test_variable, split_1, split_2);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_mci.evaluate(test_variable, split_1, split_3);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_mci.evaluate(test_variable, split_1, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_mci.evaluate(test_variable, split_2, split_3);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_mci.evaluate(test_variable, split_2, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  metric_mci.evaluate(test_variable, split_3, split_4);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), result);
  free(split_1());
  free(split_2());
  free(split_3());
  free(split_4());
}
