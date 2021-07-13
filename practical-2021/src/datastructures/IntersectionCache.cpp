#include "IntersectionCache.hpp"

    IntersectionCacheLinear::IntersectionCacheLinear(const PllPointerMap& map, const GeneralizedMetric& metric) {
      n = map.size();
      cache = std::vector<double>((n * (n + 1))/2);
      size_t k = 0;



      for(unsigned i = 0; i < n; ++i) {
          for(unsigned j = i; j < n; ++j) {
            assert(k < ((n * (n + 1))/2));
              cache[k] = metric.evaluate(i, j, map);
              ++k;
          }
      }
    }

    IntersectionCacheMatrix::IntersectionCacheMatrix(const PllPointerMap& map, const GeneralizedMetric& metric) {
      for(unsigned i = 0; i < map.size(); ++i) {
          std::vector<double> line;
          for(unsigned j = 0; j <= i; ++j) {
            line.push_back(metric.evaluate(i, j, map));
          }
          cache.push_back(line);
      }
    }
