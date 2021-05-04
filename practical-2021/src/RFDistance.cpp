#include "RFDistance.hpp"
#include "io/FileReader.hpp"

RFData RFDistance::computeRF(const std::string &data_set_path) {

  std::vector<PllSplitList> tree_splits = io::readTreeFile2(data_set_path);
  size_t tip_count = io::readTreeTipCount(data_set_path);
  size_t tree_count = tree_splits.size();

  RFData result;
  /*result.distances = std::vector<size_t>((tree_count*(tree_count-1))/2);
  //stores for every tree T_i the tree T_j which admits the smallest RF-Distance to T_i for all j < i
  std::vector<size_t> closest_tree(tree_count);
  //stores for every tree T the splits in the symmetric difference to closest_tree[T]
  std::vector<PllSplitList> D_closest;
  //stores the smallest distance encountered for the current tree so far
  size_t min_dist;
  //stores the current distance
  size_t dist;

  result.unique_count = tree_count;
  closest_tree[0] = 0;
  //D_closest.push_back(tree_splits[0]);
  for(size_t i = 1; i < tree_count; ++i){
    //D[j] stores the splits in the symmetric difference of T_j and the current tree T_i
    std::vector<PllSplitList> D;
    D.push_back(tree_splits[i].symmetricDifference(tree_splits[0]));
    //as RF distance is the cardinality of the symmetric difference
    dist = D[0].getSplitCount();
    result.distances[arrayPos(0, i, tree_count)] = dist;
    closest_tree[i] = 0;
    D_closest.push_back(PllSplitList(D[0]));
    min_dist = dist;
    if(dist == 0) --result.unique_count;
    for(size_t j = 1; j < i; ++j){
      D.push_back(D[closest_tree[j]].symmetricDifference(D_closest[j]));
      dist = D[j].getSplitCount();
      result.distances[arrayPos(j, i, tree_count)] = dist;
      if (dist < min_dist){
        closest_tree[i] = j;
        D_closest[i] = D[j];
        min_dist = dist;
        if(dist == 0) --result.unique_count;
      }
    }
  }*/


  result.unique_count = tree_splits.size();
  bool is_unique = true;
  size_t dist = 0;
  for(size_t i = 0; i < tree_splits.size(); i++){
    is_unique = true;
    for(size_t j = i+1; j < tree_splits.size(); j++){
      dist = tree_splits[i].rfDistance(tree_splits[j]);
      if (dist==0 && is_unique){
        is_unique = false;
        result.unique_count--;
      }
      result.distances.emplace_back(dist);
    }
  }

  result.relative_distances = std::vector<float>();
  for(size_t i = 0; i < result.distances.size(); ++i){
    result.relative_distances.emplace_back((float) result.distances[i] / (2*(tip_count-3)));
  }
  result.average_distance = (std::accumulate(result.distances.begin(), result.distances.end(), 0.0) / (2*(tip_count - 3))) / result.distances.size();
  result.tree_count = tree_count;
  return result;
}
