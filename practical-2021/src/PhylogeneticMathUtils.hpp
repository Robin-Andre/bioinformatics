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
  //this is an optimization which we can deploy later with mpz
  /*inline size_t truncatedDoublefactorial(size_t offset, size_t n){
    if (n <= offset) return 1;
    return n * truncatedDoublefactorial(offset, n - 2);
  }*/


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
    if ((a == 0) || (b == 0)) {
      mpq_set_ui(result, 0, 1); //TODO is 0/1 a proper result?? THis is the edge case we discussed but found no solution
      return;
    } 
    if ((a == 1) || (b == 1)) {
      mpq_set_ui(result, 1, 1); //Set result to 1/1. 
      return;
      } 
    factorialQuotient(result, ((2 * a) - 3), ((2 * b) - 3), ((2 * (a + b)) - 5));
  }
  //TODO make this private maybe?
  //TODO to get this to work with gmp we need extra tools https://github.com/linas/anant
  //right now it is a conversion to double which will cause side effects when converting really small numbers
  //MEMO aaactually since we calculate on really small numbers we could uses the inverse 
  inline double h(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    mpq_t temp_result;
    mpq_init(temp_result);
    phylogeneticProbability(temp_result, a, b);
    double conversion = mpq_get_d(temp_result);
    mpq_clear(temp_result);
    return -1 * std::log(conversion);
  }
  inline double h(const PllSplit& s) {
    return h(s.partitionSizeOf(Block_A), s.partitionSizeOf(Block_B));
  }
  //TODO this method needs a beauty session, too many ifs.
  inline void sharedPhylogeneticProbability(mpq_t result, size_t a_1, size_t b_1, size_t a_2, size_t b_2) {
    //assert(a_1 > 0 && b_1 > 0 && a_2 > 0 && b_2 > 0);
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    //in this case either identical or incompatible
    assert(a_1 != a_2); //TODO TAKE A LOOK AT THE DEFINITION AND REASSERT THAT THIS CANNOT OCCUR
    //edge cases for trivial bipartitions
    if ((a_1 == 1) || (b_1 == 1)) {
      phylogeneticProbability(result, a_2, b_2);
      return;
      }
    if ((a_2 == 1) || (b_2 == 1)) {
      phylogeneticProbability(result, a_1, b_1);
      return;
    }
    //edge case empty split
    //TODO neither of these two should not occur, debug why it did
    if ((a_1 == 0) || (b_1 == 0)) {
      phylogeneticProbability(result, a_2, b_2);
      return;
      }
    if ((a_2 == 0) || (b_2 == 0)) {
      phylogeneticProbability(result, a_1, b_1);
      return;
      }
    if(a_1 > a_2) {
      factorialQuotient(result, ((2 * a_2) - 3), ((2 * b_1) - 3), ((2 * a_1) - (2 * a_2) - 1), ((2 * (a_1 + b_1)) - 5));
      return;
    } else {
      factorialQuotient(result, ((2 * a_1) - 3), ((2 * b_2) - 3), ((2 * a_2) - (2 * a_1) - 1), ((2 * (a_2 + b_2)) - 5));
      return;
    }
  }
  //TODO either move the temporary mpq as class parameter or make sure that h is capable of dealing with mpX numbers
  inline double h(size_t a_1, size_t b_1, size_t a_2, size_t b_2){
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    mpq_t temporary_result;
    mpq_init(temporary_result);
    sharedPhylogeneticProbability(temporary_result, a_1, b_1, a_2, b_2);
    double temporary_double_holder = mpq_get_d(temporary_result);
    return -1 * std::log(temporary_double_holder);
  }
  inline double clusteringProbability(size_t count) {
    //assert(PllSplit.count > 0);
    return (1.0d * count) / PllSplit::getTipCount();
  }
  inline double clusteringProbability(const PllSplit& s, Partition block) {
      return clusteringProbability(s.partitionSizeOf(block));
  }
  inline double clusteringProbability(const PllSplit& s1, Partition block_1, const PllSplit& s2, Partition block_2) {
    //return clusteringProbability(s1.partitionSize(partition1) + s2.partitionSize(partition2) - s1.intersectionSize(s2, partition1, partition2));
    return clusteringProbability(s1.intersectionSize(s2, block_1, block_2));
  }
  inline double entropy(const PllSplit& split) {
    double p_a = clusteringProbability(split, Block_A);
    double p_b = clusteringProbability(split, Block_B);
    return -p_a * std::log(p_a) - p_b * std::log(p_b);
  }
  
}