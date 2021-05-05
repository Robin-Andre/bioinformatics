#include "RFDistance.hpp"
#include "io/FileReader.hpp"


RFData RFDistance::computeRF(const std::string &data_set_path) {
  RFData result;
  std::vector<PllSplitList> tree_splits = io::readTreeFile(data_set_path);
  result.tree_count = tree_splits.size();
  assert(result.tree_count > 0);
  result.tip_count = tree_splits[0].getSplitCount() + 3;
  assert(result.tip_count > 3);

<<<<<<< HEAD
 //std::vector<PllSplitList> tree_splits = io::readTreeFile2(data_set_path);



 std::vector<PllTree> tree_list = io::readTreeFile(data_set_path);
  std::vector<PllSplitList> tree_splits;
  tree_splits.emplace_back(PllSplitList(tree_list[0]));
  for(size_t i = 1; i < tree_list.size(); ++i) {
    PllSplitList random = PllSplitList(tree_list[i]);
    tree_splits.emplace_back(PllSplitList(random));
  }




  size_t tip_count = io::readTreeTipCount(data_set_path);
  size_t tree_count = tree_splits.size();
  /*else {
  exceptionhandeling
  }*/
=======
>>>>>>> improve_implementation

  /*result.distances = std::vector<size_t>((result.tree_count*(result.tree_count-1))/2);
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
  }*/


  result.unique_count = result.tree_count;
  bool is_unique;
  size_t dist;
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

  result.relative_distances = std::vector<float>();
  for(size_t i = 0; i < result.distances.size(); ++i){
    result.relative_distances.emplace_back((float) result.distances[i] / (2*(result.tip_count-3)));
  }
  assert(result.distances.size() == result.relative_distances.size());
  result.average_distance = std::accumulate(result.relative_distances.begin(), 
                                            result.relative_distances.end(), 0.0) / result.relative_distances.size();
  assert(result.average_distance >= 0 && result.average_distance <= 1);
  return result;
}
