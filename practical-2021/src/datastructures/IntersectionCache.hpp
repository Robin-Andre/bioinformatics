#pragma once
#include <vector>
#include "PllPointerMap.hpp"
#include "../Metric.hpp"
class IntersectionCache {
    public:

    IntersectionCache(const PllPointerMap& map, const GeneralizedMetric& metric) ;
    double access(size_t i, size_t j) const {
        return cache[pos(std::min(i, j), std::max(i, j))];
    }

    private:

  size_t pos(size_t i, size_t j)  const{
    assert(i <= j);
    size_t p = (((n*(n+1)) - ((n-i)*(n-i+1)))/2) + (j - i);
    assert(p < cache.size());
    return p;
  }


    std::vector<double> cache;
    size_t n;
};
