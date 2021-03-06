#pragma once
#include <vector>
#include "UniquePllMap.hpp"
#include "../Metric.hpp"

/*
 * Simple Class to store and access the resp. splitwise values for the metrices
 * (Matching Split Information, Shared Phylogenetic Information or
 * Mutual Clustering Information, depending on current metric)
 */

class SimilarityCache {
    public:
    /*
    * @param i The index of the first split
    * @param j The index of the second split
    * @return The value for corresponding split pair
    */
    virtual double access(size_t i, size_t j) const = 0;
    virtual ~SimilarityCache();
};

/*
 * Stores values in a linear vector, manages index calculcations for access
 */
class SimilarityCacheLinear : public SimilarityCache {
    public:
    SimilarityCacheLinear(const UniquePllMap& map, const GeneralizedMetric& metric);
    double access(size_t i, size_t j) const override {
        return cache[pos(std::min(i, j), std::max(i, j))];
    }

    virtual ~SimilarityCacheLinear() override;

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

/*
 * Stores values trivially in a two-dimensional vector
 */
class SimilarityCacheMatrix : public SimilarityCache {
    public:

    SimilarityCacheMatrix(const UniquePllMap& map, const GeneralizedMetric& metric);
    double access(size_t i, size_t j) const override{
        assert(i < cache.size() && j < cache.size());
        return cache[std::max(i, j)][std::min(i, j)];
    }

    virtual ~SimilarityCacheMatrix() override;

    private:
    std::vector<std::vector<double>> cache;
};
