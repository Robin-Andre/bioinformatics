#pragma once
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include "datastructures/PllSplits.hpp"
namespace phylomath {
    
  //TODO if gmp can do double factorial, we should use it. 
  //TODO: Is a recursive call really the way to go? Can't we improve a bit?
  inline size_t doublefactorial(size_t n) {
    if (n == 0 || n == 1) return 1; //TODO might it be smarter to check for n == 1 first, all our applications have uneven numbers
    return n * doublefactorial(n - 2);
  }
  //TODO same thing, if we get a closed form or a simple loop to unroll we'd gain performance
  inline size_t truncatedDoublefactorial(size_t offset, size_t n){
    if (n <= offset) return 1;
    return n * truncatedDoublefactorial(offset, n - 2);
  }
  //Numerators is not passed as const reference because sorting.
  inline double factorialQuotient(std::vector<size_t>& numerators, size_t denominator){
    //TODO I have no clue what this assert does :  assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    std::sort(numerators.begin(), numerators.end());

    //return doublefactorial(m) * 1.0d / truncatedDoublefactorial(M, x);
    //TODO suggested improvement
    //return (double) (doublefactorial())....etc
    return -1;
  }
  inline double factorialQuotient(size_t a, size_t b, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    size_t M = std::max(a, b);
    size_t m = std::min(a, b);
    return doublefactorial(m) * 1.0d / truncatedDoublefactorial(M, x);
  }
  //This calculates (a!!*b!!*c!!) / x!!
  inline double factorialQuotient(size_t a, size_t b, size_t c, size_t x){
    assert((a % 2 == 1) && (b % 2 == 1) && (x % 2 == 1));
    size_t M = std::max({a, b, c});
    // TODO replace with sort and maybe merge this method with the one above
    size_t m_1 = (M == a) ? ((M == b) ? c : b) : a;
    size_t m_2 = (M == b) ? c : b;
    return doublefactorial(m_1) * doublefactorial(m_2) * 1.0d / truncatedDoublefactorial(M, x); //TODO cast or multiplication
  }
  inline double phylogeneticProbability(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    if ((a == 0) || (b == 0)) return 1; //empty split TODO is this really 1 or 0? obviously a log(0) is really stupid but 0 feels right
    if ((a == 1) || (b == 1)) return 1; //trivial split
    return factorialQuotient(((2 * a) - 3), ((2 * b) - 3), ((2 * (a + b)) - 5));
  }
  //TODO make this private maybe?
  inline double h(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    return -1 * std::log(phylogeneticProbability(a, b));
  }
  //TODO don't forget the inline
  inline double h(const PllSplit& s) {
    return h(s.partitionSize(1), s.partitionSize(0));
  }

  inline double sharedPhylogeneticProbability(size_t a_1, size_t b_1, size_t a_2, size_t b_2) {
    //assert(a_1 > 0 && b_1 > 0 && a_2 > 0 && b_2 > 0);
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    //in this case either identical or incompatible
    assert(a_1 != a_2); //TODO TAKE A LOOK AT THE DEFINITION AND REASSERT THAT THIS CANNOT OCCUR
    //edge cases for trivial bipartitions
    if ((a_1 == 1) || (b_1 == 1)) return phylogeneticProbability(a_2, b_2);
    if ((a_2 == 1) || (b_2 == 1)) return phylogeneticProbability(a_1, b_1);
    //edge case empty split
    //TODO this should not occur, debug why it did
    if ((a_1 == 0) || (b_1 == 0)) return phylogeneticProbability(a_2, b_2);
    if ((a_2 == 0) || (b_2 == 0)) return phylogeneticProbability(a_1, b_1);
    if(a_1 > a_2) {
      return factorialQuotient(((2 * a_2) - 3), ((2 * b_1) - 3), ((2 * a_1) - (2 * a_2) - 1), ((2 * (a_1 + b_1)) - 5));
    } else {
      return factorialQuotient(((2 * a_1) - 3), ((2 * b_2) - 3), ((2 * a_2) - (2 * a_1) - 1), ((2 * (a_2 + b_2)) - 5));
    }
  }

  inline double h(size_t a_1, size_t b_1, size_t a_2, size_t b_2){
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    return -1 * std::log(sharedPhylogeneticProbability(a_1, b_1, a_2, b_2));
  }
  inline double clusteringProbability(size_t count) {
    //assert(PllSplit.count > 0);
    return (1.0d * count) / PllSplit::getTipCount();
  }
  inline double clusteringProbability(const PllSplit& s, partition_t partition) {
      return clusteringProbability(s.partitionSize(partition));
  }
  inline double clusteringProbability(const PllSplit& s1, partition_t partition1, const PllSplit& s2, partition_t partition2) {
    //return clusteringProbability(s1.partitionSize(partition1) + s2.partitionSize(partition2) - s1.intersectionSize(s2, partition1, partition2));
    return clusteringProbability(s1.intersectionSize(s2, partition1, partition2));
  }
  inline double entropy(const PllSplit& split) {
    double p_a = clusteringProbability(split, 1);
    double p_b = clusteringProbability(split, 0);
    return -p_a * std::log(p_a) - p_b * std::log(p_b);
  }
  
}