#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include "datastructures/PllSplits.hpp"

/**
 * Implements operations required for calculations in phylogenetics
 */

class phylomath {
private:
  static std::vector<double> ldfCache; //Note that the cache only contains log(x!!) for odd x (We dont need others)
  static std::vector<double> logCache; //Contains all logs up to Tipcount starting from 1 : log(1) log(2) ...etc.
public:

  /**
   * Fills caches for logarithms and log doublefactorials depeding
   * on the tip count set for PllSplits
   *
   */
  static void initCache() {
    initLdfCache();
    initLogCache();
  }

  /**
   * Flushes caches for logarithms and log doublefactorials
   *
   */
  static void flushCache() {
    ldfCache.clear();
    logCache.clear();
  }


  /**
   * Dual logarithm of double factorial, cache is used
   *
   * @param n
   * @return log_2(n!!)
   */
  inline static double logDoublefactorial(size_t n) {
    if (n == 0) return 0;
    if ((n - 1) / 2 < phylomath::ldfCache.size() && n % 2 == 1) {
      return phylomath::ldfCache[(n - 1) / 2];
    } else {
      return computeLogDoublefactorial(n);
    }

  }



  /**
   * Probablility that a tree on A union B contains the Split A|B
   * defined as ((2|A|−3)!!·(2|B|−3)!!) / (2(|A|+|B|)−5)!!
   *
   * @param a: |A|
   * @param b: |B|
   * @return: phylogenetic probability of A|B
   */
  inline static double phylogeneticProbability(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    assert(a > 0 && b > 0);
    if ((a == 1) || (b == 1)) {
      return 0.0;
    }
    return ldfCache[a - 2] + ldfCache[b - 2] - ldfCache[a + b - 3];
  }

  /**
   * Information content of the split Split A|B
   * defined as -1 * log_2(p), where p is the phylogenetic prob. of A|B
   *
   * @param a: |A|
   * @param b: |B|
   * @return information content of A|B
   */
  inline static double h(size_t a, size_t b) {
    assert(a + b <= PllSplit::getTipCount());
    if(a == 0 || b == 0) {
      return 0.0;
    }
    return -1 * phylogeneticProbability(a, b);
  }


  /**
   * Shared information content of splits A1|B1 and A2|B2,
   * defined as -log_2(p), where p is the phylogenetic prob. of the two splits
   * defined as (2(|A2|+ 1)−5)!!·(2(|B1|+ 1)−5)!!·(2(|A1|−|A2|+ 2)−5)!!) / (2(|A1|+|B1|)−5)!!
   * partitions are labelled such that  A2 intersect B2 =  empty set
   * @param a: |A1|
   * @param b: |B2|
   * @param
   * @return Shared information content of splits A1|B1 and A2|B2
   */

  inline static double h_shared(size_t a, size_t b) {
    size_t alltaxa = PllSplit::getTipCount();
    assert(a >= 2 && b >= 2);
    assert(a + b < alltaxa);
    size_t c;
    c = 2 * (alltaxa - a - b) - 1;
    assert(c > 0 && alltaxa > 3);
    return -1 * (ldfCache[a - 2] + ldfCache[b - 2]
              + ldfCache[alltaxa - a - b - 1] - ldfCache[alltaxa - 3]);
  }

  /**
   * Clustering probability for taxa to belong to A in a split A|B
   * defined as |A|/(|A|+|B|)    |A|+|B| = N 
   *
   * @param a: |A|
   * @return: Clustering probability for A
   */
  inline static double clusteringProbability(size_t a) {
    return static_cast<double>(a) / static_cast<double>(PllSplit::getTipCount());
  }


  /**
   * Calculates log[p_intersect / (p_1 * p_2)]. Since p_intersect = |intersect| / N; p_1 = |a1| / N; p_2 = |a2| / N 
   * this equation simplifies to log(|intersect| * N / |a_1| / |a_2|) => log(|intersect|) + log(N) - log(|a1|) - log(|a2|) 
   *
   * @param a1: |A1| the amount of taxa in split 1
   * @param a2: |A2| the amount of taxa in split 2
   * @param intersectsize: |A1 intersect A2| the amount of taxa in both split 1 and 2
   * @return: log[p_intersect / (p_1 * p_2)]. (Needed for MCI)
   */
  inline static double clusteringProbability(size_t intersectsize, size_t a1, size_t a2) {
    assert(intersectsize > 0 && a1 > 0 && a2 > 0);
    assert(intersectsize < PllSplit::getTipCount());
    assert(a1 < PllSplit::getTipCount());
    assert(a2 < PllSplit::getTipCount());
    return logCache[intersectsize - 1] + logCache[PllSplit::getTipCount() - 1]
          - logCache[a1 - 1] - logCache[a2 - 1];
  }

  /**
   * Entropy for a split A|B defined as
   * −PCl(A)*log(PCl(A))−PCl(B)*log(PCl(B))
   * where PCl is the clustering probability
   *
   * @param a: |A|
   * @param b: |B|
   * @return: Entropy for A|B
   */
  inline static double entropy(size_t a, size_t b) {
    if (a == 0 || b == 0) return 0;
    double p_a = clusteringProbability(a);
    double p_b = clusteringProbability(b);
    assert(p_a > 0 && p_b > 0);
    return -p_a * std::log2(p_a) - p_b * std::log2(p_b);
  }

private:
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
    initLogCache();
  }
  static void initLogCache() {
    logCache.clear();
    for(unsigned i = 1; i < PllSplit::getTipCount() + 1; ++i) {
      logCache.push_back(std::log2(i));
    }
  }

  inline static double computeLogDoublefactorial(size_t n) {
    if (n < 2){
      return 0;
    } else {
      return std::log2(n) + logDoublefactorial(n - 2);
    }
  }

};
