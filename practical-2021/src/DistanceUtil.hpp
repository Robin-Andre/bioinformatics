#pragma once
#include "PllSplits.hpp"
#include <tgmath.h>
#include "TriangleMatrix.hpp"
#include <vector>

enum Metric{MSI, SPI, MCI};

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

  static double h(const PllSplit& s){
    return h(s.partitionSize(1), s.partitionSize(0));
  }

  static double phylogeneticProbability(size_t a, size_t b){
    assert(a + b <= PllSplit::getTipCount());
    if ((a == 0) || (b == 0)) return 1; //empty split
    if ((a == 1) || (b == 1)) return 1; //trivial split
    return (doublefactorial((2 * a) - 3) * doublefactorial((2 * b) - 3) * 1.0d) /
      doublefactorial((2 * (a + b)) - 5);
  }

  static double MSI(const PllSplit& s1, const PllSplit& s2) {
    if (s1 ==  s2) return h(s1.partitionSize(1), s1.partitionSize(0));
    return std::max(h(s1.intersectionSize(s2, 1, 1), s1.intersectionSize(s2, 0, 0)),
                    h(s1.intersectionSize(s2, 0, 1), s1.intersectionSize(s2, 1, 0)));
  }

  static double h(size_t a_1, size_t b_1, size_t a_2, size_t b_2){
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    return -1 * std::log(sharedPhylogeneticProbability(a_1, b_1, a_2, b_2));
  }

  static double sharedPhylogeneticProbability(size_t a_1, size_t b_1, size_t a_2, size_t b_2) {
    //assert(a_1 > 0 && b_1 > 0 && a_2 > 0 && b_2 > 0);
    assert(a_1 + b_1 == PllSplit::getTipCount() && a_2 + b_2 == PllSplit::getTipCount());
    //in this case either identical or incompatible
    assert(a_1 != a_2);
    //edge cases for trivial bipartitions
    if ((a_1 == 1) || (b_1 == 1)) return phylogeneticProbability(a_2, b_2);
    if ((a_2 == 1) || (b_2 == 1)) return phylogeneticProbability(a_1, b_1);
    //edge case empty split
    if ((a_1 == 0) || (b_1 == 0)) return phylogeneticProbability(a_2, b_2);
    if ((a_2 == 0) || (b_2 == 0)) return phylogeneticProbability(a_1, b_1);
    if(a_1 > a_2) {
      return (doublefactorial((2 * a_2) - 3) * doublefactorial((2 * b_1) - 3) * doublefactorial((2 * a_1) - (2 * a_2) - 1) * 1.0d) /
        doublefactorial((2 * (a_1 + b_1)) - 5);
    } else {
      return (doublefactorial((2 * a_1) - 3) * doublefactorial((2 * b_2) - 3) * doublefactorial((2 * a_2) - (2 * a_1) - 1) * 1.0d) /
        doublefactorial((2 * (a_2 + b_2)) - 5);
    }
  }

  static double SPI(const PllSplit& s1, const PllSplit& s2){
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(s1.intersectionSize(s2, 1, 1) > 0);
    if (!s1.compatible(s2)) return 0;
    size_t a_1 = s1.partitionSize(1);
    size_t a_2 = s2.partitionSize(1);
    size_t b_1 = s1.partitionSize(0);
    size_t b_2 = s2.partitionSize(0);
    if (s1 == s2) return h(a_2, b_2);
    return h(a_1, b_1) + h(a_2, b_2) - h(a_1, b_1, a_2, b_2);

  }


  static double clusteringProbability(size_t cnt){
    assert(PllSplit::getTipCount() > 0);
    return (1.0d * cnt) / PllSplit::getTipCount();
  }

  static double clusteringProbability(const PllSplit& s, partition_t partition){
    return clusteringProbability(s.partitionSize(partition));
  }


  static double clusteringProbability(const PllSplit& s1, partition_t partition1, const PllSplit& s2, partition_t partition2) {
    //return clusteringProbability(s1.partitionSize(partition1) + s2.partitionSize(partition2) - s1.intersectionSize(s2, partition1, partition2));
    return clusteringProbability(s1.intersectionSize(s2, partition1, partition2));
  }

  static double entropy(const PllSplit& s){
    double p_a = clusteringProbability(s, 1);
    double p_b = clusteringProbability(s, 0);
    return -p_a * std::log(p_a) - p_b * std::log(p_b);
  }

  static double MCI(const PllSplit& s1, const PllSplit& s2){
    double mci = 0;
    partition_t partition1 = 0;
    do {
      partition_t partition2 = 0;
      do {
        double pcl = clusteringProbability(s1, partition1, s2, partition2);
        mci += pcl * std::log(pcl / (clusteringProbability(s1, partition1) * clusteringProbability(s2, partition2)));
        partition2 = !partition2;
      } while(partition2);
      partition1 = !partition1;
    } while(partition1);
    return mci;
  }

  static double maximumValue(const PllSplitList& first, const PllSplitList& second, Metric metric){
    assert(first.getSplits().size() == first.getSplits().size());
    size_t n = first.getSplits().size();
    double result = 0;
    switch (metric) {
      case Metric::MSI:
      case Metric::SPI:
      {
        for(size_t i = 0; i < n; ++i){
          result += h(first[i]);
          result += h(second[i]);
        }
        break;
      }
      case Metric::MCI:
      {
        for(size_t i = 0; i < n; ++i){
          result += entropy(first[i]);
          result += entropy(second[i]);
        }
        break;
      }
      default:
      {
        std::cout << "No proper metric specified" << std::endl;
        return result;
      }
    }
      return result/2;
  }

  static double distanceFromSimilarity(const PllSplitList& first, const PllSplitList& second, Metric metric, double similarity){
    return maximumValue(first, second, metric) - similarity;
  }

  static std::vector<std::vector<double>> similaritiesForSplits(const PllSplitList& first, const PllSplitList& second, Metric metric){
    assert(first.getSplits().size() == first.getSplits().size());
    size_t n = first.getSplits().size();
    std::vector<std::vector<double>>  result = std::vector<std::vector<double>>(n, std::vector<double>(n));
    for(size_t i = 0; i < n; ++i){
      for(size_t j = 0; j < n; ++j){
        switch (metric) {
          case Metric::MSI:
          {
            result[i][j] = MSI(first[i], second[j]);
            break;
          }
          case Metric::SPI:
          {
            result[i][j] = SPI(first[i], second[j]);
            assert(SPI(first[i], second[j]) != 0 || !first[i].compatible(second[j]));
            break;
          }
          case Metric::MCI:
          {
            result[i][j] = MCI(first[i], second[j]);
            break;
          }
          default:
          {
            std::cout << "No proper metric specified" << std::endl;
            return result;
          }
        }
      }
    }
    return result;
  }




};
