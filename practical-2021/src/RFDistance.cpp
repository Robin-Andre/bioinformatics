#include "RFDistance.hpp"


RFDistance::RFDistance(const std::string &data_set_path){
  std::fstream tree_file;
    tree_file.open(data_set_path, std::ios::in);
    if (tree_file.is_open()){
      std::string line;
      std::getline(tree_file, line);
      PllTree first_tree = PllTree(line);
      tip_count = first_tree.getTipCount();
      tree_splits.emplace_back(PllSplitList(first_tree));
      while(std::getline(tree_file, line)){
        PllTree tree = PllTree(line);
        tree.alignNodeIndices(first_tree);
        tree_splits.emplace_back(PllSplitList(tree));
      }
      tree_count = tree_splits.size();
      tree_file.close();
  }
}

void RFDistance::run() {
  unique_count = tree_count;
  bool is_unique = true;
  pll_split_t* s1 = new pll_split_t[tip_count-3];
  pll_split_t* s2 = new pll_split_t[tip_count-3];
  for(unsigned int i = 0; i < tree_count; i++){
    is_unique = true;
    for(unsigned int j = i+1; j < tree_count; j++){
      for(unsigned int k = 0; k < tip_count -3; k++){
        s1[k] = tree_splits[i][k].getSplit();
        s2[k] = tree_splits[j][k].getSplit();
      }
      unsigned int dist = static_cast<unsigned int(*)(pll_split_t*, pll_split_t*, unsigned int)>(&pllmod_utree_split_rf_distance)(s1, s2, tip_count);
      float normalized_dist = (float) dist / (2*(tip_count-3));
      if (dist==0 && is_unique){
        is_unique = false;
        unique_count--;
      }
      distances.emplace_back(normalized_dist);
    }
  }
  delete [] s1;
  delete [] s2;
}

std::vector<float> RFDistance::getDistances() const {
  return distances;
}

unsigned int RFDistance::getUniqueCount() const {
  return unique_count;
}

float RFDistance::getAverageDistance() const {
  return std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
}



void RFDistance::writeResults(const std::string &output_path) const {
  std::fstream distance_file;
  distance_file.open(output_path+"/distances", std::ios::out);
  if (distance_file.is_open()) {
    unsigned int k=0;
    for(unsigned int i = 0; i < tree_count; i++){
      for(unsigned int j = i+1; j < tree_count; j++){
        distance_file << i << " " << j << ": " << "0 " << distances[k] << "\n";
        k++;
        std::cout << k << std::endl;
      }
    }
    distance_file.close();
  }
}
