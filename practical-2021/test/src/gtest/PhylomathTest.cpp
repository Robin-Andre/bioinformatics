#include "gtest/gtest.h"
#include "PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"

class PhylomathTest : public testing::Test {
  protected:
    /**
     * For control of log doublefactorial
     */
    double doublefactorial(size_t n) {
      if (n < 2) {
        return 1;
      } else {
        return n * doublefactorial(n - 2);
      }
    }

    /**
     * Check correctness of log doublefactorial
     */
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


TEST_F(PhylomathTest, test_double_factorial_with_cache) {
  PllSplit::setTipCount(15);
  phylomath::initCache();
  evaluate_double_factorial(0);
  evaluate_double_factorial(1);
  evaluate_double_factorial(5);
  evaluate_double_factorial(7);
  evaluate_double_factorial(4);
  evaluate_double_factorial(17);
  evaluate_double_factorial(21);
}

TEST_F(PhylomathTest, test_double_factorial_without_cache) {
  phylomath::flushCache();
  evaluate_double_factorial(0);
  evaluate_double_factorial(1);
  evaluate_double_factorial(5);
  evaluate_double_factorial(7);
  evaluate_double_factorial(4);
  evaluate_double_factorial(17);
  evaluate_double_factorial(21);
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

/**
 * Assert that the information content of a trivial split is 0
 */
TEST_F(PhylomathTest, h_function_trivial_split) {
  PllSplit::setTipCount(4);
  PllSplit test_split = TestUtil::createSplit({0, 1, 3});
  EXPECT_EQ(test_split.h(), 0);
  free(test_split());
}

TEST_F(PhylomathTest, test_entropy) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 1, 2};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_DOUBLE_EQ(phylomath::entropy(split.partitionSizeOf(1), split.partitionSizeOf(0)), -std::log2(1.0d/2));
  free(split());
}

TEST_F(PhylomathTest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(4), 1.0d/3);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a.partitionSizeOf(1)), 5.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_a.partitionSizeOf(0)), 7.0d/12);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b.partitionSizeOf(1)), 1.0d/2);
  EXPECT_DOUBLE_EQ(phylomath::clusteringProbability(split_b.partitionSizeOf(0)), 1.0d/2);

  free(split_a());
  free(split_b());
}
