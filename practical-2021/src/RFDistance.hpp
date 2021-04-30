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
  RFDistance(const std::string &data_set_path);
  ~RFDistance(){
    for(unsigned int i=0; i < tree_splits.size(); i++){
      delete(tree_splits[i]);
    }
  }
  unsigned int getTreeCount() const {return tree_count;}
  unsigned int getTipCount() const {return tip_count;}
  void run();
  std::vector<unsigned int> getDistances() const;
  unsigned int getUniqueCount() const;
  float getAverageDistance() const;
  void writeResults(const std::string &output_path) const;

private:
  std::vector <PllSplitList*> tree_splits;
  unsigned int tree_count;
  unsigned int tip_count;
  unsigned int unique_count;
  std::vector<unsigned int> distances;

};
