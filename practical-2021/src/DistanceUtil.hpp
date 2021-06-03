#pragma once
#include "datastructures/PllSplits.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "Metric.hpp"
#include "MaximumMatcher.hpp"
#include <vector>

class DistanceUtil {


/*public:

  static double maximumValue(const PllSplitList& first, const PllSplitList& second, const Metric& metric){
    //assert(first.getSplits().size() == first.getSplits().size());
    return (metric.maximum(first, second) / 2);
  }

  static double distanceFromSimilarity(const PllSplitList& first, const PllSplitList& second, double similarity, const Metric& metric){
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

  static double distanceOf(const PllSplitList& first, const PllSplitList& second, const Metric& metric, bool normalize) {
    std::vector<std::vector<double>> similarities = similaritiesForSplits(first, second, metric);
    double similarity = MaximumMatcher::match(similarities);
    /*for (size_t k = 0; k < similarities.size(); ++k){
      for (size_t l = 0; l < similarities[k].size(); ++l){
        std::cout << similarities[k][l] << "; ";
      }
      std::cout << "|" << phylomath::h(tree_splits[i][k]);
      std::cout << std::endl;
    }
    std::cout << "SIM: " << similarity << std::endl;
    return normalize ? distanceFromSimilarity(first, second, similarity, metric) : similarity;
  }*/



};
