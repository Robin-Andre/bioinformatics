#include "RFDistance.hpp"


RFDistance::RFDistance(const std::string &data_set_path){
    std::fstream tree_file;
    tree_file.open(data_set_path, std::ios::in);
    tree_count = 0;
    if (tree_file.is_open()){
      std::string line;
      while(std::getline(tree_file, line)){
        tree_count++;
      }
      tree_file.close();
    }
    tree_file.open(data_set_path, std::ios::in);
    if (tree_file.is_open()){
      std::string line;
      std::getline(tree_file, line);
      PllTree first_tree = PllTree(line);
      tip_count = first_tree.getTipCount();
      tree_splits.reserve(tree_count);
      tree_splits[0] = new PllSplitList(first_tree);
      unsigned int i=1;
      while(std::getline(tree_file, line)){
        PllTree tree = PllTree(line);
        tree.alignNodeIndices(first_tree);
        PllSplitList* split_list = new PllSplitList(tree);
        tree_splits[i] = split_list;
        i++;
      }
      tree_file.close();
  }
}


/*void RFDistance::run() {
  HashTable<std::vector<std::pair<unsigned int, unsigned int>> H = HashTable(tree_count);
  std::map<pll_split_t, std::tuple<unsigned int, unsigned int, pll_split_t, unsigned int>> key_map;
  pll_split_t split;
  unsigned int h1;
  unsigned int h2;
  for(unsigned int i = 0; i < tree_count; i++){
    for(unsigned int j = 0; j < tip_count -3; j++){
        split = (*tree_splits[i])[j].getSplit();
        h1 = hash(split, 1);
        h2 = hash(split, 2);
        iterator = H.find(h1);
        if(H.find(h1) != H.end()){
          split2 = //get entry from iterator position
          if(key_map.count(split)){
            //collision type I
          } else if(hash(split2, 2) == h2){
            //collision type III
          } else

        }
      }
  }
}*/

/*void RFDistance::run() {
  unique_count = tree_count;
  bool is_unique = true;
  pll_split_t* s1 = new pll_split_t[tip_count-3];
  pll_split_t* s2 = new pll_split_t[tip_count-3];
  for(unsigned int i = 0; i < tree_count; i++){
    is_unique = true;
    for(unsigned int j = i+1; j < tree_count; j++){
      for(unsigned int k = 0; k < tip_count -3; k++){
        s1[k] = (*tree_splits[i])[k].getSplit();
        s2[k] = (*tree_splits[j])[k].getSplit();
      }
      unsigned int dist = pllmod_utree_split_rf_distance(s1, s2, tip_count);
      //float normalized_dist = (float) dist / (2*(tip_count-3));
      if (dist==0 && is_unique){
        is_unique = false;
        unique_count--;
      }
      distances.emplace_back(dist);
    }
  }
  delete [] s1;
  delete [] s2;
}*/


void RFDistance::run() {
  unsigned int split_count = tip_count - 3;
  unsigned int split_len   = bitv_length(tip_count);
  unique_count = tree_count;
  bool is_unique = true;
  unsigned int equal = 0;
  unsigned int s1_idx = 0,
               s2_idx = 0;
  PllSplitList* s1 = nullptr;
  PllSplitList* s2 = nullptr;
  for(unsigned int i = 0; i < tree_count; i++){
    is_unique = true;
    for(unsigned int j = i+1; j < tree_count; j++){
      s1 = tree_splits[i];
      s2 = tree_splits[j];
      equal = 0;
      s1_idx = 0;
      s2_idx = 0;
      for (s1_idx=0; s1_idx < split_count && s2_idx < split_count; s1_idx++){
        int cmp = compare_splits((*s1)[s1_idx].getSplit(), (*s2)[s2_idx].getSplit(), split_len);
        if (!cmp) {
          equal++;
          s2_idx++;
        } else {
          if (cmp > 0) {
            while(++s2_idx < split_count && (cmp = compare_splits((*s1)[s1_idx].getSplit(), (*s2)[s2_idx].getSplit(), split_len)) > 0);
            if (!cmp) {
              equal++;
            }
          }
        }
      }
      equal = 2*(tip_count - 3 - equal);
      if (equal==0 && is_unique){
        is_unique = false;
        unique_count--;
      }
      distances.emplace_back(equal);
    }
  }
}


int RFDistance::compare_splits (pll_split_t s1, pll_split_t s2, unsigned int split_len) {
  unsigned int i;
  for (i = 0; i < split_len; i++) {
    if (s1[i] != s2[i]){
      return (int) (s1[i] > s2[i]?1:-1);
    }
  }
  return 0;
}

unsigned int RFDistance::bitv_length(unsigned int bit_count) {
  unsigned int split_size = sizeof(pll_split_base_t) * 8;
  unsigned int split_offset = bit_count % split_size;
  return bit_count / split_size + (split_offset>0);
}

std::vector<unsigned int> RFDistance::getDistances() const {
  return distances;
}

unsigned int RFDistance::getUniqueCount() const {
  return unique_count;
}

float RFDistance::getAverageDistance() const {
  return (std::accumulate(distances.begin(), distances.end(), 0.0) / (2*(tip_count - 3))) / distances.size();
}



void RFDistance::writeResults(const std::string &output_path) const {
  std::fstream distance_file;
  distance_file.open(output_path+"/distances", std::ios::out);
  if (distance_file.is_open()) {
    unsigned int k=0;
    for(unsigned int i = 0; i < tree_count; i++){
      for(unsigned int j = i+1; j < tree_count; j++){
        distance_file << i << " " << j << ": " << distances[k] << "\n";
        k++;
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
