#pragma once
#include "datastructures/PllSplits.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "Metric.hpp"
#include <vector>

class DistanceUtil { //TODO remove the class and put it into a namespace (needs good suggestion) phylogenetic_

public:

  static double maximumValue(const PllSplitList& first, const PllSplitList& second, const Metric& metric) {
    //assert(first.getSplits().size() == first.getSplits().size());
    return (metric.maximum(first, second) / 2);
  }

  static double distanceFromSimilarity(const PllSplitList& first, const PllSplitList& second, const Metric& metric, double similarity){
    return maximumValue(first, second, metric) - similarity;
  }

  static std::vector<std::vector<double>> similaritiesForSplits(const PllSplitList& first, const PllSplitList& second, const Metric& metric){
    assert(first.getSplits().size() == first.getSplits().size());
    size_t n = first.getSplits().size();
    std::vector<std::vector<double>>  result = std::vector<std::vector<double>>(n, std::vector<double>(n));
    for(size_t i = 0; i < n; ++i){
      for(size_t j = 0; j < n; ++j){
        result[i][j] = metric.evaluate(first[i], second[j]);
      }
    }
    return result;
  }




};
