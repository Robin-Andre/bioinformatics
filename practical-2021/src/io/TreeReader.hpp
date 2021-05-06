#pragma once
//#include<vector>
#include <fstream>
//#include <iostream>
#include "../PllTree.hpp"
//#include "../RFDistance.hpp"
//#include "CommandLineOptions.hpp"
class TreeReader {

public:

static std::vector<PllSplitList> readTreeFile(const std::string& filepath) {
  std::vector<PllSplitList> pll_list;
  std::ifstream file(filepath);
  if(file.is_open()) {
    std::string line;
    std::getline(file, line);
    PllTree first_tree = PllTree(line);
    PllSplit::setSplitLen(PllSplit::computeSplitLen(first_tree.getTipCount()));
    pll_list.emplace_back(PllSplitList(first_tree));
    while(std::getline(file, line)) {
      PllTree tree_from_line = PllTree(line);
      assert(first_tree.getTipCount() == tree_from_line.getTipCount());
      tree_from_line.alignNodeIndices(first_tree);
      pll_list.emplace_back(PllSplitList(tree_from_line));
    }
  }
  file.close();
  return pll_list;
}
};
