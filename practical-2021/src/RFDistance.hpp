#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllTree.hpp"
#include "RFData.hpp"
#include "io/TreeReader.hpp"
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <numeric>
#include <bitset>
#include <algorithm>
#include <cfloat>


class RFDistance {
public:
  static RFData computeRF(const std::vector<PllTree>& trees) {
    size_t tree_count = trees.size();
    assert(tree_count > 0);
    size_t tip_count = trees[0].getTipCount();
    assert(tip_count > 3);
    PllSplit::setTipCount(tip_count);
    std::vector<PllSplitList> tree_splits;
    for(PllTree tree :  trees){
      tree_splits.emplace_back(PllSplitList(tree));
    }

    size_t dist = 0;
    size_t unique_count = tree_count;
    bool is_unique = true;
    std::vector<double> distances;
    for(size_t i = 0; i < tree_count; i++){
      is_unique = true;
      for(size_t j = i+1; j < tree_count; j++){
        dist = tree_splits[i].rfDistance(tree_splits[j]);
        assert(dist <= 2*(tip_count - 3));
        if (dist==0 && is_unique){
          is_unique = false;
          unique_count--;
        }
        distances.emplace_back(dist);
      }
    }
    return RFData(tree_count, tip_count, unique_count, distances, true);
  }
};
