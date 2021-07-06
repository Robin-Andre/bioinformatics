#pragma once
#include <vector>
#include "PllPointerMap.hpp"
#include "../Metric.hpp"
class IntersectionCache {
    public:

    IntersectionCache(const PllPointerMap& map, const GeneralizedMetric& metric) ;
    double access(size_t i, size_t j) const {
        return cache[i][j];
    }

    private:

 
    std::vector<std::vector<double>> cache;
};