#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../PllTree.hpp"
namespace io {
static std::vector<PllTree> readTreeFile(const std::string& filepath) {
  std::vector<PllTree> tree_vector;
  std::ifstream file(filepath);
  if(file.is_open()) {
      std::string line;
      while(std::getline(file, line)) {
          PllTree test(line);
          std::cout << line << "\n";
          tree_vector.push_back(test);
      }
  }
  file.close();
  //tree_vector.push_back(PllTree("roflcopter"));
  return tree_vector;

}
static size_t getDistanceFromString(const std::string &line) {
    std::istringstream iss (line);
    std::string item;
    size_t i = 0;
    while (std::getline(iss, item, ' ') && i < 2) {
        i++;

    }
    return std::stoi(item);
}
static void writeResultFile(const std::string& filepath){

}
static std::vector<size_t> readDistances(std::string data_set_name) {
    std::vector<size_t> distances;
    std::fstream res_file;
    res_file.open("../test/res/reference_results/" + data_set_name + "/RAxML_RF-Distances.0"  ,std::ios::in);
    if (res_file.is_open()){
      std::string line;
      std::vector<std::string> parts;
      while(std::getline(res_file, line)){
        distances.push_back(getDistanceFromString(line));
      }
      res_file.close(); //close the file object.
  }
  return distances;
}


}//namespace io