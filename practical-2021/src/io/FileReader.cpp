#include "FileReader.hpp"
#include <fstream>
#include <iostream>

std::vector<PllTree>& readTreeFile(const std::string& filepath) {
  std::vector<PllTree> tree_vector;
  std::ifstream file(filepath);
  if(file.is_open()) {
      std::string line;
      while(std::getline(file, line)) {
          std::cout << line << "\n";
      }
  }
  //tree_vector.push_back(PllTree("roflcopter"));
  return tree_vector;
}
static void writeResultFile(const std::string& filepath) {
    //TODO build and fix 
}

