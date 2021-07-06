#include "IntersectionCache.hpp"

    IntersectionCache::IntersectionCache(const PllPointerMap& map, const GeneralizedMetric& metric) {
      cache = std::vector(map.size(), std::vector<double>(map.size()));
      for(unsigned i = 0; i < cache.size(); ++i) {
          //Theoretically we could make this a triangular matrix and make the access well formed i.e. A < B
          for(unsigned j = 0; j < cache[i].size(); ++j) {
              cache[i][j] = metric.evaluate(i, j, map);            
          }
      }
    }
