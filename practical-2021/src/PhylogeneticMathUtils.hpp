#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include "datastructures/PllSplits.hpp"

class phylomath {
private:
  static std::vector<double> ldfCache; //Note that the cache only contains log(x!!) for odd x (We dont need others)
public:

  inline static double computeLogDoublefactorial(size_t n) {
    if (n < 2){
      return 0;
    } else {
      return std::log2(n) + logDoublefactorial(n - 2);
    }
  }
  static void initLdfCache() {
    size_t upper_bound = (PllSplit::getTipCount() * 2) - 1;
    size_t lower_bound;
    if(ldfCache.size() == 0){
      ldfCache.push_back(0);
      lower_bound = 3;
    } else {
      lower_bound = (ldfCache.size() * 2) - 1;
    }
    double prev = ldfCache.back();
    for (size_t i = lower_bound; i <= upper_bound; i+=2){
      prev = prev + std::log2(i);
      phylomath::ldfCache.push_back(prev);
    }
  }

  static void flushLdfCache() {
    ldfCache.clear();
  }


  //Check for cache hit should be removed, make sure, that cache always hits
  //@Robin: I will try to move all logDF calls up one method, if the methods are fine no
  // cache miss should ever occur
  inline static double logDoublefactorial(size_t n) {
    if (n == 0) return 0;
    if ((n - 1) / 2 < phylomath::ldfCache.size() && n % 2 == 1) {
      return phylomath::ldfCache[(n - 1) / 2];
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

    return ldfCache[a - 2] + ldfCache[b - 2] - ldfCache[a + b - 3];
    //return logFactorialQuotient(((2 * a) - 3), ((2 * b) - 3), ((2 * (a + b)) - 5));
  }


  inline static double h(size_t a, size_t b) {
    assert(a + b <= PllSplit::getTipCount());
    if(a == 0 || b == 0) {
      return 0.0d;
    }
    return -1 * phylogeneticProbability(a, b);
  }

  //This method is a mockup of the calculation of phylogenetic probability of two splits
  //Requires the size of the two partitions (A or B) which need to be compatible in the first place
  //The calculation will work even if they are not compatible but the result is entirely useless
  inline static double h(size_t taxa_partition1, size_t taxa_partition2, size_t alltaxa) {
    assert(taxa_partition1 >= 2 && taxa_partition2 >= 2);
    assert(taxa_partition1 + taxa_partition2 < alltaxa);
    size_t c;
    c = 2 * (alltaxa - taxa_partition1 - taxa_partition2) - 1;
    assert(c > 0 && alltaxa > 3);
    return -1 * (ldfCache[taxa_partition1 - 2] + ldfCache[taxa_partition2 - 2]
              + ldfCache[alltaxa - taxa_partition1 - taxa_partition2 - 1] - ldfCache[alltaxa - 3]);
    //return -1 * logFactorialQuotient(a, b, c, x);

  }

  inline static double clusteringProbability(size_t count) {
    return static_cast<double>(count) / static_cast<double>(PllSplit::getTipCount());
  }

  inline static double clusteringProbability(const PllSplit* s1, Partition block_1, const PllSplit* s2, Partition block_2) {
    return clusteringProbability(s1->intersectionSize(*s2, block_1, block_2));
  }
  inline static double entropy(size_t a, size_t b) {
    if (a == 0 || b == 0) return 0;
    double p_a = clusteringProbability(a);
    double p_b = clusteringProbability(b);
    assert(p_a != 0 && p_b != 0);
    return -p_a * std::log2(p_a) - p_b * std::log2(p_b);
  }

};
