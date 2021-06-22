#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include "datastructures/PllSplits.hpp"
class phylomath {
private:
  static std::vector<double> ldfCache;

public:

  inline static double computeLogDoublefactorial(size_t n) {
    if (n < 2){
      return 0;
    } else {
      return std::log2(n) + computeLogDoublefactorial(n - 2);
    }
  }

  static void initLdfCache() {
    for (size_t i = ldfCache.size(); i <= PllSplit::getTipCount(); ++i){
      phylomath::ldfCache.push_back(computeLogDoublefactorial(i));
    }
  }

  /* This is the "apparent" GMP double factorial function I dislike the style of void functions changing the
     value of their parameters but in this case it is inevitable due to the implementation of gmp.*/

  inline static double logDoublefactorial(size_t n) {
    if (n < phylomath::ldfCache.size()){
      return phylomath::ldfCache[n];
    } else {
      return computeLogDoublefactorial(n);
    }

  }



  //This calculates log((a!!*b!!*c!!) / x!!)
  //TODO this method needs a beauty session
  inline static double logFactorialQuotient(size_t a, size_t b, size_t c, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    return logDoublefactorial(a) + logDoublefactorial(b) + logDoublefactorial(c) - logDoublefactorial(x);
  }

  inline static double logFactorialQuotient(size_t a, size_t b, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    return logFactorialQuotient(a, b, 1, x);
  }



  inline static double phylogeneticProbability(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    assert(a > 0 && b > 0);
    if ((a == 1) || (b == 1)) {
      return 0.0d;
    }
    return logFactorialQuotient(((2 * a) - 3), ((2 * b) - 3), ((2 * (a + b)) - 5));
  }


  inline static double h(size_t a, size_t b) {
    assert(a + b <= PllSplit::getTipCount());
    if(a == 0 || b == 0) {
      return 0.0d;
    }
    return -1 * phylogeneticProbability(a, b);;
  }

  inline static double h(const PllSplit& s) {
    //There should never be a partition where one block is empty
    assert(s.partitionSizeOf(Block_A) > 0 && s.partitionSizeOf(Block_B) > 0);
    return h(s.partitionSizeOf(Block_A), s.partitionSizeOf(Block_B));
  }
  //This method is a mockup of the calculation of phylogenetic probability of two splits
  //Requires the size of the two partitions (A or B) which need to be compatible in the first place
  //The calculation will work even if they are not compatible but the result is entirely useless
  inline static double h(size_t taxa_partition1, size_t taxa_partition2, size_t alltaxa) {
    assert(taxa_partition1 >= 2 && taxa_partition2 >= 2);
    assert(taxa_partition1 + taxa_partition2 < alltaxa);
    size_t a, b, c, x;
    a = 2 * taxa_partition1 - 3;
    b = 2 * taxa_partition2 - 3;
    c = 2 * (alltaxa - taxa_partition1 - taxa_partition2) - 1;
    x = 2 * alltaxa - 5;
    assert(a > 0 && b > 0 && c > 0 && x > 0);
    return -1 * logFactorialQuotient(a, b, c, x);;

  }

  inline static double clusteringProbability(size_t count) {
    return static_cast<double>(count) / static_cast<double>(PllSplit::getTipCount());
  }
  inline static double clusteringProbability(const PllSplit& s, Partition block) {
      return clusteringProbability(s.partitionSizeOf(block));
  }
  inline static double clusteringProbability(const PllSplit& s1, Partition block_1, const PllSplit& s2, Partition block_2) {
    return clusteringProbability(s1.intersectionSize(s2, block_1, block_2));
  }
  inline static double entropy(const PllSplit& split) {
    double p_a = clusteringProbability(split, Block_A);
    double p_b = clusteringProbability(split, Block_B);
    assert(p_a != 0 && p_b != 0);
    return -p_a * std::log2(p_a) - p_b * std::log2(p_b);
  }

};
