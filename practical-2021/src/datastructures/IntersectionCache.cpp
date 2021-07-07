#include "IntersectionCache.hpp"

    IntersectionCache::IntersectionCache(const PllPointerMap& map, const GeneralizedMetric& metric) {
      n = map.size();
      cache = std::vector<double>((n * (n + 1))/2);
      size_t k = 0;
      for(unsigned i = 0; i < cache.size(); ++i) {
          for(unsigned j = i; j < n; ++j) {
              cache[k] = metric.evaluate(i, j, map);
              ++k;
          }
      }
    }
