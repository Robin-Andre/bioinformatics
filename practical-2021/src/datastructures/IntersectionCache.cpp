#include "IntersectionCache.hpp"

    IntersectionCache::IntersectionCache(const PllPointerMap& map, const GeneralizedMetric& metric) {
      cache = std::vector(map.size(), std::vector<double>(map.size()));
      for(unsigned i = 0; i < cache.size(); ++i) {
          for(unsigned j = i + 1; j < cache[i].size(); ++j) {
              cache[i][j] = metric.evaluate(i, j, map);
              cache[j][i] = cache[i][j];            
          }
      }
      for(unsigned i = 0; i < cache.size(); ++i) {

        cache[i][i] = metric.evaluate(i, i, map);           
          
      }
    }
