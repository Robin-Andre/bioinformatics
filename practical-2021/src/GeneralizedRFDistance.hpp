#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllTree.hpp"
#include "io/IOData.hpp"
#include "Metric.hpp"
#include <vector>
#include <iostream>

//TODO without attributes it might not need to be a class
class GeneralizedRFDistance {
public:

  static io::IOData computeDistances(const std::vector<PllTree>& trees, const Metric& metric, Mode mode) {
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);
    PllSplit::setTipCount(trees[0].getTipCount());
    phylomath::initLdfCache();
    std::vector<PllSplitList> tree_splits;
    for(PllTree tree :  trees){
      assert(tree.getTipCount() == PllSplit::getTipCount());
      tree_splits.emplace_back(PllSplitList(tree));
    }
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
        double dist = metric.distanceOf(tree_splits[i], tree_splits[j], mode);
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

};
