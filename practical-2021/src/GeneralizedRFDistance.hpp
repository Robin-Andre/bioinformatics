#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "DistanceUtil.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"
#include "datastructures/TriangleMatrix.hpp"
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
    TriangleMatrix<double> distances = TriangleMatrix<double>(tree_count, false);
    double similarity = 0;
    double dist = 0;
    size_t unique_count = trees.size();
    bool is_unique = true;
    for(size_t i = 0; i < trees.size(); ++i){
      is_unique = true;
      for(size_t j = i + 1; j < trees.size(); ++j){
        std::vector<std::vector<double>> similarities = DistanceUtil::similaritiesForSplits(tree_splits[i], tree_splits[j], metric);
        for (size_t k = 0; k < similarities.size(); ++k){
          for (size_t l = 0; l < similarities[k].size(); ++l){
            //std::cout << similarities[k][l] << "; ";
          }
          //std::cout << "|" << phylomath::h(tree_splits[i][k]);
          //std::cout << std::endl;
        }
        similarity = MaximumMatcher::match(similarities);
        //std::cout << "SIM: " << similarity << std::endl;

        dist = normalize ? DistanceUtil::distanceFromSimilarity(tree_splits[i], tree_splits[j], metric, similarity) : similarity;
        if (dist == 0 && is_unique){
          is_unique = false;
          --unique_count;
        }
        distances.set(i, j, dist);
      }
    }
    return RFData(tree_count, tip_count, unique_count, distances.getAsVector(), false);
  }

};
