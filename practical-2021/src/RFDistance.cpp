#include "RFDistance.hpp"


RFData RFDistance::computeRF(const std::vector<PllTree>& trees) {
  RFData result;
  result.tree_count = trees.size();
  assert(result.tree_count > 0);
  result.tip_count = trees[0].getTipCount();
  assert(result.tip_count > 3);
  PllSplit::setTipCount(result.tip_count);
  std::vector<PllSplitList> tree_splits;
  for(PllTree tree :  trees){
    tree_splits.emplace_back(PllSplitList(tree));
  }


  result.distances = std::vector<double>((result.tree_count*(result.tree_count-1))/2);
  //stores for every tree T_i the tree T_j which admits the smallest RF-Distance to T_i for all j < i
  std::vector<size_t> closest_tree(result.tree_count);
  //stores for every tree T the splits in the symmetric difference to closest_tree[T]
  std::vector<PllSplitList> D_closest;
  //stores the smallest distance encountered for the current tree so far
  size_t min_dist;
  //stores the current distance
  size_t dist;

  result.unique_count = result.tree_count;
  closest_tree[0] = 0;
  D_closest.push_back(tree_splits[0]);
  for(size_t i = 1; i < result.tree_count; ++i){
    //D[j] stores the splits in the symmetric difference of T_j and the current tree T_i
    std::vector<PllSplitList> D;
    D.push_back(tree_splits[i].symmetricDifference(tree_splits[0]));
    //as RF distance is the cardinality of the symmetric difference
    dist = D[0].getSplitCount();
    assert(dist >= 0 && dist <= 2*(result.tip_count - 3));
    result.distances[arrayPos(0, i, result.tree_count)] = dist;
    closest_tree[i] = 0;
    D_closest.push_back(PllSplitList(D[0]));
    assert(D_closest.size() == i+1);
    min_dist = dist;
    if(dist == 0) --result.unique_count;
    for(size_t j = 1; j < i; ++j){
      D.push_back(D[closest_tree[j]].symmetricDifference(D_closest[j]));
      assert(D.size() == j+1);
      dist = D[j].getSplitCount();
      assert(dist >= 0 && dist <= 2*(result.tip_count - 3));
      result.distances[arrayPos(j, i, result.tree_count)] = dist;
      if (dist < min_dist){
        closest_tree[i] = j;
        PllSplitList D2 = D[j];
        D_closest[i] = D[j];
        min_dist = dist;
        if(dist == 0) --result.unique_count;
      }
    }
  }


  /*result.unique_count = result.tree_count;
  bool is_unique = true;
  size_t dist = 0;
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
  assert(result.distances.size() == (result.tree_count * (result.tree_count - 1)) / 2);*/


  result.average_distance = (std::accumulate(result.distances.begin(),
                                            result.distances.end(), 0.0) / (2*(result.tip_count-3))) / result.distances.size();
  assert(result.average_distance >= 0 && result.average_distance <= 1);
  return result;
}
