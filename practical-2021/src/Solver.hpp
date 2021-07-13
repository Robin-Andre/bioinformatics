#pragma once
#include <vector>
#include "datastructures/PllSplits.hpp"
class Solver {

  public:

    static double distanceFromSimilarity( const PllSplitList& first,
                                          const PllSplitList& second,
                                          double similarity,
                                          Mode mode,
                                          const GeneralizedMetric& metric) {
      double max_value = metric.maximum(first, second);
      assert(max_value != 0);
      assert((max_value - 2 * similarity) > -0.000001);
      double dist = std::max(0.0, max_value - 2 * similarity);
      return (mode == RELATIVE) ? dist : (dist / max_value);
    }

    static void similaritiesForSplits(const PllSplitList& first, const PllSplitList& second, const PllPointerMap& map,
         std::vector<std::vector<double>>* result, const IntersectionCache& cache){
      assert(first.getSplits().size() == first.getSplits().size());
      size_t n = first.getSplits().size();
      assert(result->size() == n);
      for(size_t i = 0; i < n; ++i){
        assert((*result)[i].size() == n);
        for(size_t j = 0; j < n; ++j) {
          (*result)[i][j] = cache.access(first[i], second[j]);
        }
      }
    }

}; //class solver
