#pragma once
#include "PllSplits.hpp"
#include <tgmath.h>

class DistanceUtil {

public:

  static size_t doublefactorial(size_t n) {
    if (n == 0 || n == 1) return 1;
    return n * doublefactorial(n - 2);
  }


  static double h(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    return -1 * std::log(phylogeneticProbability(a, b));
  }

  static double phylogeneticProbability(size_t a, size_t b){
    return (doublefactorial((2 * a) - 3) * doublefactorial((2 * b) - 3) * 1.0d) /
      doublefactorial((2 * (a + b)) - 5);
  }

  static double MSI(PllSplit s1, PllSplit s2) {
    return std::max(h(s1.intersectionSize(s2, false, false), s1.intersectionSize(s2, true, true)),
                    h(s1.intersectionSize(s2, true, false), s1.intersectionSize(s2, false, true)));
  }

  static double h(size_t a_1, size_t b_1, size_t a_2, size_t b_2){
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    return -1 * std::log(sharedPhylogeneticProbability(a_1, b_1, a_2, b_2));
  }

  static double sharedPhylogeneticProbability(size_t a_1, size_t b_1, size_t a_2, size_t b_2) {
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    if(a_1 > a_2) {
      return (doublefactorial((2 * a_2) - 3) * doublefactorial((2 * b_1) - 3) * doublefactorial((2 * a_1) - (2 * a_2) - 1) * 1.0d) /
        doublefactorial((2 * (a_1 + b_1)) - 5);
    } else {
      return (doublefactorial((2 * a_1) - 3) * doublefactorial((2 * b_2) - 3) * doublefactorial((2 * a_2) - (2 * a_1) - 1) * 1.0d) /
        doublefactorial((2 * (a_2 + b_2)) - 5);
    }
  }

  static double SPI(PllSplit s1, PllSplit s2){
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(s1.intersectionSize(s2, false, false) > 0);
    if (!s1.compatible(s2)) return 0;
    size_t a_1 = s1.popcount();
    size_t a_2 = s2.popcount();
    size_t b_1 = PllSplit::getTipCount() - a_1;
    size_t b_2 = PllSplit::getTipCount() - a_2;
    return h(a_1, b_1) + h(a_2, b_2) - h(a_1, b_1, a_2, b_2);

  }



  static double clusteringProbability(PllSplit s, bool invert){
    size_t a = s.popcount();
    if (!invert) return a / PllSplit::getTipCount();
    else return 1 - (a / PllSplit::getTipCount());
  }

  static double clusteringInformationContent(PllSplit s, bool invert) {
    return -1 * std::log(clusteringProbability(s, invert));
  }

  static double entropy(PllSplit s){
    double p_a = clusteringProbability(s, false);
    double p_b = clusteringProbability(s, true);
    return -p_a * std::log(p_a) - p_b * std::log(p_b);
  }

  static double MCI(PllSplit s1, PllSplit s2){
    //?!
  }






};
