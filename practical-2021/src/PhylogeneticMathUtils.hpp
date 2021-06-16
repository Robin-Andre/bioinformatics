#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include "datastructures/PllSplits.hpp"
namespace phylomath {

  /* This is the "apparent" GMP double factorial function I dislike the style of void functions changing the
     value of their parameters but in this case it is inevitable due to the implementation of gmp.*/

  inline void logdoublefactorial(mpfr_t acc, size_t n, size_t offset) {
    if (n <= offset){
      mpfr_set_ui(acc, 0, RND);
    } else {
      mpfr_t s;
      mpfr_init_set_ui(s, n, RND);
      mpfr_log2(s, s, RND);
      logdoublefactorial(acc, n-2, offset);
      mpfr_add(acc, acc, s, RND);
      mpfr_clear(s);
    }
  }



  inline void doublefactorial(mpz_t result, size_t n) {
    mpz_2fac_ui(result, n);
  }

  //This calculates (a!!*b!!*c!!) / x!!
  //TODO this method needs a beauty session
  inline void factorialQuotient(mpq_t result, size_t a, size_t b, size_t c, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    /*TODO optimization, move this allocation out of this method, maybe make it a class
      because setting a mpz to 0 could be better than to reinitialize every single time
      a factorialQuotient is calculated. (Which might be often)
    */
    mpz_t fac_a, fac_b, fac_c, fac_x;
    mpz_init(fac_a);
    mpz_init(fac_b);
    mpz_init(fac_c);
    mpz_init(fac_x);
    doublefactorial(fac_a, a);
    doublefactorial(fac_b, b);
    doublefactorial(fac_c, c);
    doublefactorial(fac_x, x);
    mpz_mul(fac_a, fac_a, fac_b);
    mpz_mul(fac_a, fac_a, fac_c);
    mpq_set_num(result, fac_a);
    mpq_set_den(result, fac_x);
    mpq_canonicalize(result); //TODO could be a slowdown, needs evaluation

    //TODO optimization, own class and clear methods in the destructor to save calls
    mpz_clear(fac_a);
    mpz_clear(fac_b);
    mpz_clear(fac_c);
    mpz_clear(fac_x);
  }
  inline void factorialQuotient(mpq_t result, size_t a, size_t b, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    factorialQuotient(result, a, b, 1, x);
  }




  static void logFactorialQuotient(mpfr_t result, size_t a, size_t b, size_t c, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    size_t M, m_1, m_2;
    if (a >= b && a >= c){
      M = a;
      m_1 = b;
      m_2 = c;
    } else if (b >= a && b >= c){
      M = b;
      m_1 = a;
      m_2 = c;
    } else {
      assert(c >= a && c >= b);
      M = c;
      m_1 = a;
      m_2 = b;
    }

    mpfr_t counter_1, counter_2, denominator;
    mpfr_init_set_ui(counter_1, 0, RND);
    mpfr_init_set_ui(counter_2, 0, RND);
    mpfr_init_set_ui(denominator, 0, RND);
    logdoublefactorial(counter_1, m_1, 1);
    logdoublefactorial(counter_2, m_2, 1);
    logdoublefactorial(denominator, x, M);
    mpfr_add(counter_1, counter_1, counter_2, RND);
    mpfr_sub(result, counter_1, denominator, RND);
    mpfr_clear(counter_1);
    mpfr_clear(counter_2);
    mpfr_clear(denominator);
  }

  static void logFactorialQuotient(mpfr_t result, size_t a, size_t b, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    logFactorialQuotient(result, a, b, 1, x);
  }







  inline void phylogeneticProbability(mpfr_t result, size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    assert(a > 0 && b > 0);
    if ((a == 1) || (b == 1)) {
      mpfr_set_ui(result, 0, RND); //Set result to 1/1.
      return;
      }
    logFactorialQuotient(result, ((2 * a) - 3), ((2 * b) - 3), ((2 * (a + b)) - 5));
    mpfr_mul_si(result, result, -1, RND);


  }
  //TODO to get this to work with gmp we need extra tools https://github.com/linas/anant
  //right now it is a conversion to double which will cause side effects when converting really small numbers
  //MEMO actually since we calculate on really small numbers we could theoretically use the inverse
  inline void h(mpfr_t result, size_t a, size_t b) {
    assert(a + b <= PllSplit::getTipCount());
    if(a == 0 || b == 0) {
      mpfr_set_ui(result, 0, RND);
      return;
    }
    phylogeneticProbability(result, a, b);

  }
  inline void h(mpfr_t result, const PllSplit& s) {
    //There should never be a partition where one block is empty
    assert(s.partitionSizeOf(Block_A) > 0 && s.partitionSizeOf(Block_B) > 0);
    h(result, s.partitionSizeOf(Block_A), s.partitionSizeOf(Block_B));
  }
  //This method is a mockup of the calculation of phylogenetic probability of two splits
  //Requires the size of the two partitions (A or B) which need to be compatible in the first place
  //The calculation will work even if they are not compatible but the result is entirely useless
  inline void h(mpfr_t result, size_t taxa_partition1, size_t taxa_partition2, size_t alltaxa) {
    assert(taxa_partition1 >= 2 && taxa_partition2 >= 2);
    assert(taxa_partition1 + taxa_partition2 < alltaxa);
    /* If the partitions are compatible and the splits nonequal then
    there has to be at least one taxa which is in neither partition */
    size_t a, b, c, x;
    a = 2 * taxa_partition1 - 3;
    b = 2 * taxa_partition2 - 3;
    //(c) calculates the remaining tree (which should NOT be empty)
    c = 2 * (alltaxa - taxa_partition1 - taxa_partition2) - 1;

    x = 2 * alltaxa - 5;
    assert(a > 0 && b > 0 && c > 0 && x > 0);
    logFactorialQuotient(result, a, b, c, x);
    mpfr_mul_si(result, result, -1.0, RND);

  }

  inline void clusteringProbability(mpfr_t result, size_t count) {
    mpfr_set_ui(result, count, RND);
    mpfr_div_ui(result, result, PllSplit::getTipCount(), RND);
  }
  inline void clusteringProbability(mpfr_t result, const PllSplit& s, Partition block) {
      clusteringProbability(result, s.partitionSizeOf(block));
  }
  inline void clusteringProbability(mpfr_t result, const PllSplit& s1, Partition block_1, const PllSplit& s2, Partition block_2) {
    clusteringProbability(result, s1.intersectionSize(s2, block_1, block_2));
  }
  inline void entropy(mpfr_t result, const PllSplit& split) {
    mpfr_t p_a, p_b, p_a_log, p_b_log;
    mpfr_init_set_ui(p_a, 0, RND);
    mpfr_init_set_ui(p_b, 0, RND);
    mpfr_init_set_ui(p_a_log, 0, RND);
    mpfr_init_set_ui(p_b_log, 0, RND);
    clusteringProbability(p_a, split, Block_A);
    clusteringProbability(p_b, split, Block_B);
    mpfr_log2(p_a_log, p_a, RND);
    mpfr_log2(p_b_log, p_b, RND);
    mpfr_mul(p_a, p_a, p_a_log, RND);
    mpfr_mul(p_b, p_b, p_b_log, RND);
    mpfr_mul_si(p_a, p_a, -1, RND);
    mpfr_sub(result, p_a, p_b, RND);
    mpfr_clear(p_a);
    mpfr_clear(p_b);
    mpfr_clear(p_a_log);
    mpfr_clear(p_b_log);
  }

}
