#include "RFDistance.hpp"


RFData RFDistance::computeRF(const std::string &data_set_path) {

  std::vector<PllSplitList> tree_splits;
  size_t tip_count;
  std::fstream tree_file;
  tree_file.open(data_set_path, std::ios::in);
  size_t tree_count = 0;
  if (tree_file.is_open()){
    std::string line;
    while(std::getline(tree_file, line)){
      ++tree_count;
    }
    tree_file.close();
  }
  tree_file.open(data_set_path, std::ios::in);
  if (tree_file.is_open()){
    std::string line;
    std::getline(tree_file, line);
    PllTree first_tree = PllTree(line);
    tip_count = first_tree.getTipCount();
    PllSplitList first_split = PllSplitList(first_tree);
    tree_splits.emplace_back(first_split);
    size_t i=1;
    while(std::getline(tree_file, line)){
      PllTree tree = PllTree(line);
      tree.alignNodeIndices(first_tree);
      PllSplitList split_list = PllSplitList(tree);
      tree_splits.emplace_back(split_list);
      ++i;
    }
    tree_file.close();
  }
  /*else {
  exceptionhandeling
  }*/

  RFData result;
  result.distances = std::vector<size_t>((tree_count*(tree_count-1))/2);
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
  D_closest.push_back(tree_splits[0]);
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
  }

  result.relative_distances = std::vector<float>();
  for(size_t i = 0; i < result.distances.size(); ++i){
    result.relative_distances.emplace_back(result.distances[i] / (2*(tip_count-3)));
  }
  result.average_distance = (std::accumulate(result.distances.begin(), result.distances.end(), 0.0) / (2*(tip_count - 3))) / result.distances.size();

  return result;
}
