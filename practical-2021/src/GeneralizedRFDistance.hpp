#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "RFData.hpp"
#include "Metric.hpp"
#include <vector>
#include <iostream>

//TODO without attributes it might not need to be a class
class GeneralizedRFDistance {
public:

  static RFData computeDistances(const std::vector<PllTree>& trees, const Metric& metric, bool normalize) {
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);
    PllSplit::setTipCount(trees[0].getTipCount());
    std::vector<PllSplitList> tree_splits;
    for(PllTree tree :  trees){
      assert(tree.getTipCount() == PllSplit::getTipCount());
      tree_splits.emplace_back(PllSplitList(tree));
    }
    std::vector<double> distances;
    double similarity = 0;
    double dist = 0;
    size_t unique_count = trees.size();
    bool is_unique = true;
    for(size_t i = 0; i < trees.size(); ++i){
      is_unique = true;
      for(size_t j = i + 1; j < trees.size(); ++j){
        dist = metric.distanceOf(tree_splits[i], tree_splits[j], normalize);
        if (dist == 0 && is_unique){
          is_unique = false;
          --unique_count;
        }
        distances.emplace_back(dist);
      }
    }
    return RFData(tree_count, tip_count, unique_count, distances, metric.name() == "RF");
  }

};
