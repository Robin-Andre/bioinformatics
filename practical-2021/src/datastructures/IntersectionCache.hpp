#pragma once
#include <vector>
#include "PllPointerMap.hpp"
class Metric;
class IntersectionCache {
    public:

    IntersectionCache(const PllPointerMap& map, const Metric& metric) ;
    double access(size_t i, size_t j) const {
        return cache[i][j];
    }

    private:

 
    std::vector<std::vector<double>> cache;
};