#pragma once

#include <fstream>
#include "../datastructures/PllTree.hpp"

/*
 * Class to read several trees in newick format from file
 */
class TreeReader {

public:

  /*
   * @param filepath: path of file with newick strings
   * @return vector with PllTrees read from specified file
   */
static std::vector<PllTree> readTreeFile(const std::string& filepath) {
  std::vector<PllTree> trees;
  std::ifstream file(filepath);
  if(file.is_open()) {
    std::string line;
    std::getline(file, line);
    PllTree first_tree = PllTree(line);
    trees.emplace_back(first_tree);
    while(std::getline(file, line)) {
      PllTree tree_from_line = PllTree(line, first_tree);
      assert(first_tree.getTipCount() == tree_from_line.getTipCount());
      trees.emplace_back(tree_from_line);
    }
  }
  file.close();
  return trees;
}
};
