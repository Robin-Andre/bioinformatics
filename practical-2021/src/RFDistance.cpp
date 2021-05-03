#include "RFDistance.hpp"


void RFDistance::run(const std::string &data_set_path) {
  std::vector<PllSplitList> tree_splits;
  std::fstream tree_file;
  tree_file.open(data_set_path, std::ios::in);
  tree_count = 0;
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
  distances = std::vector<size_t>((tree_count*(tree_count-1))/2);
  //stores for every tree T_i the tree T_j which admits the smallest RF-Distance to T_i for all j < i
  std::vector<size_t> closest_tree(tree_count);
  //stores for every tree T the splits in the symmetric difference to closest_tree[T]
  std::vector<PllSplitList> D_closest;
  //stores the smallest distance encountered for the current tree so far
  size_t min_dist;
  //stores the current distance
  size_t dist;

  unique_count = tree_count;
  closest_tree[0] = 0;
  D_closest.push_back(tree_splits[0]);
  for(size_t i = 1; i < tree_count; ++i){
    //D[j] stores the splits in the symmetric difference of T_j and the current tree T_i
    std::vector<PllSplitList> D;
    D.push_back(tree_splits[i].symmetricDifference(tree_splits[0]));
    //as RF distance is the cardinality of the symmetric difference
    dist = D[0].getSplitCount();
    distances[getPos(0, i)] = dist;
    closest_tree[i] = 0;
    D_closest.push_back(PllSplitList(D[0]));
    min_dist = dist;
    if(dist == 0) --unique_count;
    for(size_t j = 1; j < i; ++j){
      D.push_back(D[closest_tree[j]].symmetricDifference(D_closest[j]));
      dist = D[j].getSplitCount();
      distances[getPos(j, i)] = dist;
      if (dist < min_dist){
        closest_tree[i] = j;
        D_closest[i] = D[j];
        min_dist = dist;
        if(dist == 0) --unique_count;
      }
    }
  }
}





void RFDistance::writeResults(const std::string &output_path) const {
  std::fstream distance_file;
  distance_file.open(output_path+"/distances", std::ios::out);
  if (distance_file.is_open()) {
    size_t k=0;
    for(size_t i = 0; i < tree_count; ++i){
      for(size_t j = i+1; j < tree_count; ++j){
        distance_file << i << " " << j << ": " << distances[k] << "\n";
        ++k;
      }
    }
    distance_file.close();
  }

  std::fstream info_file;
  info_file.open(output_path+"/info", std::ios::out);
  if (info_file.is_open()) {
    info_file << "Found " << tree_count << " trees\n";
    info_file << "Average relative RF in this set: " << getAverageDistance() << "\n";
    info_file << "Number of unique trees in this tree set: " << unique_count << "\n";
    info_file.close();
  }

}
