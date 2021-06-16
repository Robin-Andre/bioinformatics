#include "gtest/gtest.h"
#include "../../../src/PhylogeneticMathUtils.hpp"
#include "../TestUtil.hpp"
#include <gmp.h>
#include <mpfr.h>

class PhylomathTest : public testing::Test {
  protected:
  mpfr_t test_variable;
  mpfr_t result_variable;
  mpz_t test_variable_int;
  mpz_t result_variable_int;
  mpq_t test_variable_rational;
  mpq_t result_variable_rational;
  virtual void SetUp() {
    mpz_init(test_variable_int);
    mpz_init(result_variable_int);
    mpq_init(test_variable_rational);
    mpq_init(result_variable_rational);
    mpfr_init_set_ui(test_variable, 0, RND);
    mpfr_init_set_ui(result_variable, 0, RND);
  }
  virtual void TearDown() {
    mpz_clear(test_variable_int);
    mpz_clear(result_variable_int);
    mpq_clear(test_variable_rational);
    mpq_clear(result_variable_rational);
    mpfr_clear(test_variable);
    mpfr_clear(result_variable);
  }
  //Evaluates the double factorial and compares it to a string of base 10
  void evaluate_double_factorial(size_t x, const std::string& result_string) {
    phylomath::doublefactorial(test_variable_int, x);
    //gmp_printf("test: %Zd \n", test_variable);
    mpz_set_str(result_variable_int, result_string.c_str(), 10);
    //gmp_printf("result: %Zd \n", result_variable);
    EXPECT_EQ(mpz_cmp(test_variable_int, result_variable_int), 0);
  }
  void evaluate_phylogenetic_probability(size_t a, size_t b, const std::string& fraction) {
    phylomath::phylogeneticProbability(test_variable, a, b);
    mpfr_set_str(result_variable, fraction.c_str(), 10, RND); //Base 10


  }

  void evaluate_log_double_factorial(size_t x) {
    phylomath::logdoublefactorial(test_variable, x, 1);
    double ld = mpfr_get_d(test_variable, RND);
    mpz_2fac_ui(test_variable_int, x);
    double d = mpz_get_d(test_variable_int);
    EXPECT_DOUBLE_EQ(ld, std::log2(d));
  }

  void evaluate_log_factorial_quotient(size_t a, size_t b, size_t c, size_t x) {
    phylomath::logFactorialQuotient(test_variable, a, b, c, x);
    double lq = mpfr_get_d(test_variable, RND);

    phylomath::factorialQuotient(test_variable_rational, a, b, c, x);
    double q = mpq_get_d(test_variable_rational);
    EXPECT_DOUBLE_EQ(lq, std::log2(q));
  }

};



TEST_F(PhylomathTest, test_double_factorial) {
  evaluate_double_factorial(0, "1");
  evaluate_double_factorial(1, "1");
  evaluate_double_factorial(2, "2");
  evaluate_double_factorial(3, "3");
  evaluate_double_factorial(4, "8");
  evaluate_double_factorial(6, "48");
  evaluate_double_factorial(41, "13113070457687988603440625");
  evaluate_double_factorial(60, "284813089515958324736640819941867520000000");
  evaluate_double_factorial(65, "7297912393562140321551086320493608726062890625");
}
//TODO I am not entirely sure what I should test, especially since mpqs are really finnicy and
//double factorial is tested also this test will fail if normalization is remeoved in factorial quotient
TEST_F(PhylomathTest, test_factorial_quotient) {
  phylomath::factorialQuotient(test_variable_rational, 1001, 1, 1003);
  mpq_set_str(result_variable_rational, "1/1003", 10);
  int comp = mpq_cmp(test_variable_rational, result_variable_rational);
  EXPECT_EQ(comp, 0);
}
//9 <- (this eats 2nominators) * 7 * 5 * 3 <- this eats the last so result should be 1/35
TEST_F(PhylomathTest, test_factorial_quotient2) {
  mpq_t quotientresult;
  mpq_init(quotientresult);
  phylomath::factorialQuotient(quotientresult, 3, 3, 3, 9);
  mpz_t thirtyfive;
  mpz_init(thirtyfive);
  mpz_set_ui(thirtyfive, 35);
  int comp = mpz_cmp(thirtyfive, mpq_denref(quotientresult));
  EXPECT_EQ(comp, 0);
}



TEST_F(PhylomathTest, test_phylogenetic_probability) {

  PllSplit::setTipCount(8);
  evaluate_phylogenetic_probability(2, 3, "1/5");
  evaluate_phylogenetic_probability(2, 4, "1/7");
  evaluate_phylogenetic_probability(3, 4, "1/21");
  evaluate_phylogenetic_probability(2, 2, "1/3");
  evaluate_phylogenetic_probability(3, 3, "3/35");
  evaluate_phylogenetic_probability(4, 4, "5/231");
  PllSplit::setTipCount(24);
  evaluate_phylogenetic_probability(2,22, "1/43");
}
//The probability of a trivial split is 1, the value of h should be 0 as in -log(1) == 0
TEST_F(PhylomathTest, h_function_trivial_split) {
  PllSplit::setTipCount(4);
  PllSplit test_split = TestUtil::createSplit({0, 1, 3});
  phylomath::h(test_variable, test_split);
  EXPECT_EQ(mpfr_get_d(test_variable, RND), 0);
  free(test_split());
}

TEST_F(PhylomathTest, test_entropy) {
  PllSplit::setTipCount(6);
  std::vector<size_t> part1 = {0, 1, 2};
  PllSplit split = TestUtil::createSplit(part1);
  phylomath::entropy(test_variable, split);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), -std::log2(1.0d/2));
  free(split());
}

TEST_F(PhylomathTest, test_clustering_probability) {
  PllSplit::setTipCount(12);
  phylomath::clusteringProbability(test_variable, 4);
  mpfr_set_ui(result_variable, 1, RND);
  mpfr_div_ui(result_variable, result_variable, 3, RND);
  EXPECT_EQ(mpfr_cmp(result_variable, test_variable), 0);
  std::vector<size_t> part1_a = {0, 1, 2, 3, 4};
  PllSplit split_a = TestUtil::createSplit(part1_a);
  std::vector<size_t> part1_b = {0, 3, 4, 5, 6, 7};
  PllSplit split_b = TestUtil::createSplit(part1_b);
  phylomath::clusteringProbability(test_variable, split_a, Block_A);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 5.0d/12);
  phylomath::clusteringProbability(test_variable, split_a, Block_B);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 7.0d/12);
  phylomath::clusteringProbability(test_variable, split_b, Block_A);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0d/2);
  phylomath::clusteringProbability(test_variable, split_b, Block_B);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0d/2);

  phylomath::clusteringProbability(test_variable, split_a, Block_A, split_b, Block_A);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0d/4);
  phylomath::clusteringProbability(test_variable, split_a, Block_A, split_b, Block_B);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0d/6);
  phylomath::clusteringProbability(test_variable, split_a, Block_B, split_b, Block_A);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0d/4);
  phylomath::clusteringProbability(test_variable, split_a, Block_B, split_b, Block_B);
  EXPECT_DOUBLE_EQ(mpfr_get_d(test_variable, RND), 1.0d/3);

  free(split_a());
  free(split_b());
}

TEST_F(PhylomathTest, test_log_double_factorial) {
  evaluate_log_double_factorial(0);
  evaluate_log_double_factorial(1);
  evaluate_log_double_factorial(2);
  evaluate_log_double_factorial(3);
  evaluate_log_double_factorial(4);
  evaluate_log_double_factorial(5);
  evaluate_log_double_factorial(6);
  evaluate_log_double_factorial(21);
}

TEST_F(PhylomathTest, test_log_factorial_quotient) {
  evaluate_log_factorial_quotient(3, 3, 3, 9);
  evaluate_log_factorial_quotient(3, 5, 7, 17);
}
