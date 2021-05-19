#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "PllSplits.hpp"
#include "DistanceUtil.hpp"
#include "MaximumMatcher.hpp"
#include "TriangleMatrix.hpp"
#include <vector>
#include <iostream>


class GeneralizedRFDistance {
public:

  TriangleMatrix<double> computeDistances(const std::vector<PllTree>& trees, Metric metric) {
    PllSplit::setTipCount(trees[0].getTipCount());
    std::vector<PllSplitList> tree_splits;
    for(PllTree tree :  trees){
      assert(tree.getTipCount() == PllSplit::getTipCount());
      tree_splits.emplace_back(PllSplitList(tree));
    }
    TriangleMatrix<double> result = TriangleMatrix<double>(trees.size(), false);
    double similarity = 0;
    double dist = 0;
    size_t unique_count = trees.size();
    bool is_unique = true;
    for(size_t i = 0; i < trees.size(); i++){
      is_unique = true;
      for(size_t j = i+1; j < trees.size(); j++){
        std::vector<std::vector<double>> similarities = DistanceUtil::similaritiesForSplits(tree_splits[i], tree_splits[j], metric);
        /*for (size_t k = 0; k < similarities.size(); ++k){
          for (size_t l = 0; l < similarities[k].size(); ++l){
            std::cout << similarities[k][l] << "; ";
          }
          std::cout << std::endl;
        }*/
        similarity = MaximumMatcher::match(similarities);
        //std::cout << similarity << std::endl;
        dist = DistanceUtil::distanceFromSimilarity(tree_splits[i], tree_splits[j], metric, similarity);
        if (dist == 0 && is_unique){
          is_unique = false;
          unique_count--;
        }
        result.set(i, j, dist);
      }
    }
    return result;
  }
private:

};
