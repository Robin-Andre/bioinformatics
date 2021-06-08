#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include <gmp.h>
#include "datastructures/PllSplits.hpp"
namespace phylomath {

  /* This is the "apparent" GMP double factorial function I dislike the style of void functions changing the
     value of their parameters but in this case it is inevitable due to the implementation of gmp.*/

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







  inline void phylogeneticProbability(mpq_t result, size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    assert(a > 0 && b > 0);
    if ((a == 1) || (b == 1)) {
      mpq_set_ui(result, 1, 1); //Set result to 1/1.
      return;
      }
    factorialQuotient(result, ((2 * a) - 3), ((2 * b) - 3), ((2 * (a + b)) - 5));
  }
  //TODO to get this to work with gmp we need extra tools https://github.com/linas/anant
  //right now it is a conversion to double which will cause side effects when converting really small numbers
  //MEMO actually since we calculate on really small numbers we could theoretically use the inverse
  inline double h(size_t a, size_t b) {
    assert(a + b <= PllSplit::getTipCount());
    if(a == 0 || b == 0) {
      return 0;
    }
    mpq_t temp_result;
    mpq_init(temp_result);
    phylogeneticProbability(temp_result, a, b);
    double conversion = mpq_get_d(temp_result);
    mpq_clear(temp_result);
    return -1 * std::log2(conversion);
  }
  inline double h(const PllSplit& s) {
    //There should never be a partition where one block is empty
    assert(s.partitionSizeOf(Block_A) > 0 && s.partitionSizeOf(Block_B) > 0);
    return h(s.partitionSizeOf(Block_A), s.partitionSizeOf(Block_B));
  }
  //This method is a mockup of the calculation of phylogenetic probability of two splits
  //Requires the size of the two partitions (A or B) which need to be compatible in the first place
  //The calculation will work even if they are not compatible but the result is entirely useless
  inline double h(size_t taxa_partition1, size_t taxa_partition2, size_t alltaxa) {
    assert(taxa_partition1 >= 2 && taxa_partition2 >= 2);
    assert(taxa_partition1 + taxa_partition2 < alltaxa);
    /* If the partitions are compatible and the splits nonequal then
    there has to be at least one taxa which is in neither partition */
    mpq_t temporary_result;
    mpq_init(temporary_result);
    size_t a, b, c, x;
    a = 2 * taxa_partition1 - 3;
    b = 2 * taxa_partition2 - 3;
    //(c) calculates the remaining tree (which should NOT be empty)
    c = 2 * (alltaxa - taxa_partition1 - taxa_partition2) - 1;

    x = 2 * alltaxa - 5;
    assert(a > 0 && b > 0 && c > 0 && x > 0);
    factorialQuotient(temporary_result, a, b, c, x);
    double temporary_double_holder = mpq_get_d(temporary_result);
    mpq_clear(temporary_result);
    return -1 * std::log2(temporary_double_holder);

  }

  inline double clusteringProbability(size_t count) {
    return static_cast<double>(count) / static_cast<double>(PllSplit::getTipCount());
  }
  inline double clusteringProbability(const PllSplit& s, Partition block) {
      return clusteringProbability(s.partitionSizeOf(block));
  }
  inline double clusteringProbability(const PllSplit& s1, Partition block_1, const PllSplit& s2, Partition block_2) {
    return clusteringProbability(s1.intersectionSize(s2, block_1, block_2));
  }
  inline double entropy(const PllSplit& split) {
    double p_a = clusteringProbability(split, Block_A);
    double p_b = clusteringProbability(split, Block_B);
    assert(p_a != 0 && p_b != 0);
    return -p_a * std::log2(p_a) - p_b * std::log2(p_b);
  }

}
