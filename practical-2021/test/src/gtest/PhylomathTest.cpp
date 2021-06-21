#include "gtest/gtest.h"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"

class PhylomathTest : public testing::Test {
  protected:

    double doublefactorial(size_t n) {
      if (n < 2) {
        return 1;
      } else {
        return n * doublefactorial(n - 2);
      }
    }

    void evaluate_double_factorial(size_t n){
      EXPECT_DOUBLE_EQ(std::log2(doublefactorial(n)), phylomath::logDoublefactorial(n));
    }
  void evaluate_phylogenetic_probability(size_t a, size_t b, double fraction) {
    EXPECT_NEAR(phylomath::phylogeneticProbability(a, b), std::log2(fraction), 0.000000000001);
  }

};



TEST_F(PhylomathTest, test_double_factorial) {
  evaluate_double_factorial(0);
  evaluate_double_factorial(1);
  evaluate_double_factorial(2);
  evaluate_double_factorial(3);
  evaluate_double_factorial(4);
  evaluate_double_factorial(6);
  evaluate_double_factorial(20);
}

TEST_F(PhylomathTest, test_factorial_quotient) {
  EXPECT_NEAR(phylomath::logFactorialQuotient(1001, 1, 1003), std::log2(1.0d/1003), 0.000000000001);
}
TEST_F(PhylomathTest, test_factorial_quotient2) {
  EXPECT_DOUBLE_EQ(phylomath::logFactorialQuotient(3, 3, 3, 9), std::log2(1.0d/35));
}



TEST_F(PhylomathTest, test_phylogenetic_probability) {

  PllSplit::setTipCount(8);
  evaluate_phylogenetic_probability(2, 3, 1.0d/5);
  evaluate_phylogenetic_probability(2, 4, 1.0d/7);
  evaluate_phylogenetic_probability(3, 4, 1.0d/21);
  evaluate_phylogenetic_probability(2, 2, 1.0d/3);
  evaluate_phylogenetic_probability(3, 3, 3.0d/35);
  evaluate_phylogenetic_probability(4, 4, 5.0d/231);
  PllSplit::setTipCount(24);
  evaluate_phylogenetic_probability(2,22, 1.0d/43);
}
//The probability of a trivial split is 1, the value of h should be 0 as in -log(1) == 0
TEST_F(PhylomathTest, h_function_trivial_split) {
  PllSplit::setTipCount(4);
  PllSplit test_split = TestUtil::createSplit({0, 1, 3});
  EXPECT_EQ(phylomath::h(test_split), 0);
  free(test_split());
}

TEST_F(PhylomathTest, test_entropy) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 1, 2};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_DOUBLE_EQ(phylomath::entropy(split), -std::log2(1.0d/2));
  free(split());
}

TEST_F(PhylomathTest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(4), 1.0d/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_A), 5.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_B), 7.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b, Block_A), 1.0d/2);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b, Block_B), 1.0d/2);

  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_A, split_b, Block_A), 1.0d/4);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_A, split_b, Block_B), 1.0d/6);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_B, split_b, Block_A), 1.0d/4);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a, Block_B, split_b, Block_B), 1.0d/3);

  free(split_a());
  free(split_b());
}
