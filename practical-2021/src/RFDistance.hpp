#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "PllSplits.hpp"
#include "PllTree.hpp"
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <numeric>
#include <bitset>


class RFDistance {
public:
  RFDistance(){};
  ~RFDistance(){};
  void run(const std::string &data_set_path);

  size_t getTreeCount() const {return tree_count;}
  size_t getTipCount() const {return tip_count;}

  std::vector<size_t> getDistances() const {return distances;}
  std::vector<float> getRelativeDistances() const {
    std::vector<float> relative_distances;
    for(size_t i = 0; i < distances.size(); ++i){
      relative_distances.emplace_back(distances[i] / (2*(getTipCount()-3)));
    }
    return relative_distances;
  }
  size_t getUniqueCount() const {return unique_count;}
  float getAverageDistance() const {return (std::accumulate(distances.begin(), distances.end(), 0.0) / (2*(tip_count - 3))) / distances.size();}
  void writeResults(const std::string &output_path) const;

private:
  //std::vector <PllSplitList*> tree_splits;
  size_t tree_count;
  size_t tip_count;
  size_t unique_count;
  std::vector<size_t> distances;

  //determines position in linear array, make sure that i < j
  size_t getPos(size_t i, size_t j) const {
    size_t offset =(i*(2*tree_count-i-1))/2;
    return offset + (j - i - 1);
  }

};
