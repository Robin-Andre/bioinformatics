#pragma once
#include "PllSplits.hpp"
#include "RFDistance.hpp"
#include "PllTree.hpp"
#include <string>
#include <iostream>

int main() {
  std::vector<std::string> tree_strings{
      "(a, b, (c, d));",
      "(a, d, (c, b));",
  };
  std::vector<PllTree> tree_list;
  for (auto t : tree_strings) { tree_list.emplace_back(t); }

  std::vector<PllSplitList> splits_list;
  std::cout << " HELLO WORLD" << std::endl;
  for (auto &t : tree_list) {
    t.alignNodeIndices(*tree_list.begin());
    splits_list.emplace_back(t);
  }
}
