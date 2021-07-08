#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllTree.hpp"
#include "datastructures/PllPointerMap.hpp"
#include "datastructures/IntersectionCache.hpp"
#include "io/IOData.hpp"
#include "Metric.hpp"
#include "Solver.hpp"
#include <vector>
#include <iostream>
#include <queue>
#include <future>
class GeneralizedRFDistance {
public:
  static io::IOData computeGeneralizedDistances(const std::vector<PllTree>& trees, const GeneralizedMetric& metric, Mode mode) {
    PllSplit::setTipCount(trees[0].getTipCount());
    //test.initializePointerSplitLists(trees);
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);

    phylomath::initLdfCache();
    PllPointerMap _map(trees);
    IntersectionCacheMatrix pairwise_cache(_map, metric);


    std::vector<PllSplitList>& tree_splits = _map.vectors();

    io::IOData result;
    result.mode = ModeString[mode];
    result.metric = metric.name();
    result.mean_dst = 0;
    result.number_of_unique_trees = tree_count;
    result.pairwise_distance_mtx = std::vector<std::vector<double>>(tree_count, std::vector<double>());
    size_t dist_count = 0;
    size_t split_count = tree_splits[0].getSplits().size();
    std::vector<std::vector<double>> similarities = std::vector<std::vector<double>>(split_count, std::vector<double>(split_count));
    std::vector<std::future<double>> futures;
    for(size_t i = 0; i < tree_count; ++i){
      futures.clear();
      for(size_t j = i; j < tree_count; ++j){
        Solver::similaritiesForSplits(tree_splits[i], tree_splits[j], _map, &similarities, pairwise_cache);
        futures.push_back(std::async(MaximumMatcher::match, similarities));
      }
      bool is_unique = true;
      for(size_t j = i; j < tree_count; ++j){
        double dist = futures[j-i].get();
        assert(dist >= 0.0);
        dist = (mode == SIMILARITY) ? dist : Solver::distanceFromSimilarity(tree_splits[i], tree_splits[j], dist, mode, metric);
        //TODO: Check near 0 because of numerical issues
        if (i != j && dist == 0 && is_unique){
          is_unique = false;
          --result.number_of_unique_trees;
        }
        result.pairwise_distance_mtx[j].emplace_back(dist);
        if(i != j){
          result.mean_dst += dist;
          ++dist_count;
        }
      }
    }
    result.mean_dst = (dist_count == 0) ? 0 : result.mean_dst  / static_cast<double>(dist_count);
    assert(result.number_of_unique_trees >= 1);
    return result;
  }


  static io::IOData computeRFDistances(const std::vector<PllTree>& trees, const RFMetric& metric, Mode mode) {
    PllSplit::setTipCount(trees[0].getTipCount());
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);
    PllPointerMap map(trees);
    std::vector<PllSplitList>& tree_splits = map.vectors();
    io::IOData result;
    result.mode = ModeString[mode];
    result.metric = metric.name();
    result.mean_dst = 0;
    result.number_of_unique_trees = tree_count;
    result.pairwise_distance_mtx = std::vector<std::vector<double>>(tree_count, std::vector<double>());
    size_t dist_count = 0;
    for(size_t i = 0; i < tree_count; ++i){
      bool is_unique = true;
      for(size_t j = i; j < tree_count; ++j){
        double dist = metric.distanceOfRF(tree_splits[i], tree_splits[j], mode, map);
        assert(dist >= 0.0);
        //TODO: Check near 0 because of numerical issues
        if (i != j && dist == 0 && is_unique){
          is_unique = false;
          --result.number_of_unique_trees;
        }
        result.pairwise_distance_mtx[j].emplace_back(dist);
        if(i != j){
          result.mean_dst += dist;
          ++dist_count;
        }
      }
    }
    result.mean_dst = (dist_count == 0) ? 0 : result.mean_dst  / static_cast<double>(dist_count);
    assert(result.number_of_unique_trees >= 1);
    return result;
  }

  private:
  PllPointerMap _map;


};
