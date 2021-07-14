#pragma once
#include <vector>
#include "PllPointerMap.hpp"
#include "../Metric.hpp"

class IntersectionCache {
    public:
    virtual double access(size_t i, size_t j) const = 0;
    virtual ~IntersectionCache() {}
};


class IntersectionCacheLinear : public IntersectionCache {
    public:
    IntersectionCacheLinear(const PllPointerMap& map, const GeneralizedMetric& metric);
    double access(size_t i, size_t j) const override {
        return cache[pos(std::min(i, j), std::max(i, j))];
    }

    virtual ~IntersectionCacheLinear() override {}

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

class IntersectionCacheMatrix : public IntersectionCache {
    public:

    IntersectionCacheMatrix(const PllPointerMap& map, const GeneralizedMetric& metric);
    double access(size_t i, size_t j) const override{
        assert(i < cache.size() && j < cache.size());
        return cache[std::max(i, j)][std::min(i, j)];
    }

    virtual ~IntersectionCacheMatrix() override {}

    private:
    std::vector<std::vector<double>> cache;
};
