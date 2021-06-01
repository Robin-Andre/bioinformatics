#include "RFDistance.hpp"


RFData RFDistance::computeRF(const std::vector<PllTree>& trees) {
  RFData result(trees.size(), trees[0].getTipCount());
  assert(result.tree_count > 0);
  assert(result.tip_count > 3);
  PllSplit::setTipCount(result.tip_count);
  std::vector<PllSplitList> tree_splits;
  for(PllTree tree :  trees){
    tree_splits.emplace_back(PllSplitList(tree));
  }

  size_t dist = 0;
  result.unique_count = result.tree_count;
  bool is_unique = true;
  for(size_t i = 0; i < result.tree_count; i++){
    is_unique = true;
    for(size_t j = i+1; j < result.tree_count; j++){
      dist = tree_splits[i].rfDistance(tree_splits[j]);
      assert(dist <= 2*(result.tip_count - 3));
      if (dist==0 && is_unique){
        is_unique = false;
        result.unique_count--;
      }
      result.distances.emplace_back(dist);
    }
  }
  assert(result.distances.size() == (result.tree_count * (result.tree_count - 1)) / 2);
  result.average_distance = (std::accumulate(result.distances.begin(),
                                            result.distances.end(), 0.0) / (2*(result.tip_count-3))) / result.distances.size();
  assert(result.average_distance >= 0 && result.average_distance <= 1);
  return result;
}
