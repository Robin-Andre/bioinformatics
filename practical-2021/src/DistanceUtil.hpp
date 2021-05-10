#pragma once
#include "PllSplits.hpp"
#include <tgmath.h>

class DistanceUtil {

public:

  static size_t doublefactorial(size_t n) {
    if (n == 0 || n == 1) return 1;
    return n * doublefactorial(n - 2);
  }

  /*static double phylogeneticProbability(const PllSplit& s){
    size_t a = s.popCount();
    size_t b = tip_count - a;
    return (doublefactorial((2 * a) - 3) * doublefactorial((2 * b) - 3)) /
      doublefactorial((2 * tip_count) - 5);
  }*/

  static double informationContent(size_t a, size_t b){
    return -1 * std::log(phylogeneticProbability(a, b));
  }

  static double phylogeneticProbability(size_t a, size_t b){
    return (doublefactorial((2 * a) - 3) * doublefactorial((2 * b) - 3) * 1.0d) /
      doublefactorial((2 * (a + b)) - 5);
  }

  static double MSI(PllSplit s1, PllSplit s2) {
    return std::max(informationContent(s1.intersectcount(s2, false, false), s1.intersectcount(s2, true, true)),
                    informationContent(s1.intersectcount(s2, true, false), s1.intersectcount(s2, false, true)));
  }

  //static double sharedPhylogeneticProbability(PllSplit s1, PllSplit s2)



  /*static double phylogeneticProbability(PllSplit s1, PllSplit s2){
    size_t a1 = s1.popCount();
    size_t b1 = tip_count - a1;
    size_t a2 = s2.popCount();
    size_t b2 = tip_count - a2;
    //Here it must be determined, which parts of the splits overlap!!!!!!!!!!!!!!!!!!
    return (doublefactorial((2 * a2) - 3) *
            doublefactorial((2 * b1) - 3) *
            doublefactorial((2 * (a1 - a2)) - 1)) /
            doublefactorial((2 * tip_count) - 5);
  }

  static double informationContent(PllSplit s1, PllSplit s2){
    return -1 * std::log(phylogeneticProbability(s1, s2));
  }

  static double SPI(PllSplit s1, PllSplit s2){
    if(! s1.compatibel(s2)) return 0;
    return informationContent(s1) + informationContent(s2) - informationContent(s1, s2);
  }

  static double clusteringProbability(PllSplit s, bool partition){
    size_t a = s.popcount();
    if (partition) return a / tip_count;
    else return 1 - (a / tip_count);
  }

  static double clusteringInformationContent(PllSplit s, bool partition) {
    return -1 * std::log(clusteringProbability(s, partition));
  }

  static double entropy(PllSplit s){
    double p_a = clusteringProbability(s, 1);
    double p_b = clusteringProbability(s, 0);
    return -p_a*std::log(p_a) - p_b*std::log(p_b);
  }

  static double MCI(PllSplit s1, PllSplit s2){
    std::vector<PllSplit> splits = {s1, s2}
    //?!
  }*/






};
