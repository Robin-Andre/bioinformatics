#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../PllTree.hpp"
#include "../RFDistance.hpp"
#include "CommandLineOptions.hpp"
namespace io {
static std::vector<PllTree> readTreeFile(const std::string& filepath) {
  std::vector<PllTree> tree_vector;
  std::ifstream file(filepath);
  if(file.is_open()) {
      std::string line;
      std::getline(file, line);
      PllTree align_tree(line);
      tree_vector.emplace_back(align_tree);
      while(std::getline(file, line)) {
          PllTree another_tree = PllTree(line);
          another_tree.alignNodeIndices(align_tree);
          tree_vector.emplace_back(another_tree);
      }
  }
  file.close();
  //tree_vector.push_back(PllTree("roflcopter"));
  return tree_vector;
}

static void writeOutput(const RFData& result, const Config& config) {
  if(config.output_file_path.size() <= 1) {
    std::cerr << "The Output was not properly specified...exiting";
    exit(1);
  }
  std::ofstream out_file(config.output_file_path);
  unsigned k = 0;
  for(unsigned i = 0; i < result.tree_count; ++i) {

    for(unsigned j = i + 1; j < result.tree_count; ++j) {
      out_file << i << " " << j << " " << result.distances[k] << " " << result.relative_distances[k] << std::endl;
      ++k;
    }
  }
  out_file.close();
  //for(const auto &e : result.distances) out_file << e << "\n";
  std::cout << "Result File written at: "<< config.output_file_path << std::endl;
  std::cout << "Unique count: " << result.unique_count
            << " Average distance: " << result.average_distance
            << " Tree count: " << result.tree_count << "\n";
  /*for(size_t i = 0; i < result.distances.size(); ++i) {
    std::cout << result.distances[i] << "\n";
  }*/
}

static std::vector<PllSplitList> readTreeFile(const std::string& filepath) {
  std::vector<PllSplitList> pll_list;
  std::ifstream file(filepath);
  if(file.is_open()) {
    std::string line;
    std::getline(file, line);
    PllTree first_tree = PllTree(line);
    pll_list.emplace_back(PllSplitList(first_tree));
    while(std::getline(file, line)) {
      PllTree tree_from_line = PllTree(line);
      tree_from_line.alignNodeIndices(first_tree);
      pll_list.emplace_back(PllSplitList(tree_from_line));
    }
  }
  file.close();
  return pll_list;
}


static int getDistanceFromString(const std::string &line) {
    std::istringstream iss (line);
    std::string item;
    size_t i = 0;
    while (std::getline(iss, item, ' ') && i < 2) {
        i++;

    }
    return std::stoi(item);
}

static std::vector<size_t> readDistances(const std::string& data_set_name) {
    std::vector<size_t> distances;
    std::fstream res_file;
    res_file.open("../test/res/reference_results/" + data_set_name + "/RAxML_RF-Distances.0"  ,std::ios::in);
    if (res_file.is_open()){
      std::string line;
      while(std::getline(res_file, line)){
        distances.push_back((size_t) getDistanceFromString(line));
      }
      res_file.close(); //close the file object.
  }
  return distances;
}
static const std::string readFromInfoFile(const std::string& data_set_name,const std::string& prefix) {
   std::string result;
   std::fstream res_file;
   res_file.open("../test/res/reference_results/" + data_set_name + "/RAxML_info.0"  ,std::ios::in);
   if (res_file.is_open()){
     std::string line;
     while(std::getline(res_file, line)){
       auto match = std::mismatch(prefix.begin(), prefix.end(), line.begin());
       if (match.first == prefix.end())
       {
         std::istringstream iss (line);
         std::istringstream iss2 (prefix);
         std::string dummy;
         while (std::getline(iss2, dummy, ' ')){
           std::getline(iss, result, ' ');
         }
         std::getline(iss, result, ' ');
         break;
       }
     }
     res_file.close(); //close the file object.
 }
 return result;
}
/*static size_t readTreeCount(const std::string& data_set_name) {
  return std::stoi(readFromInfoFile(data_set_name, "Found "));
}*/

static int readUniqueTreeCount(const std::string& data_set_name) {
  return std::stoi(readFromInfoFile(data_set_name, "Number of unique trees in this tree set: "));
}


static float readAverageDistance(const std::string& data_set_name) {
  return std::stof(readFromInfoFile(data_set_name, "Average relative RF in this set: "));
}

/*static RFData readRAxML(std::string data_set_name) {
  RFData raxml_result;
  return raxml_result;
}*/







}//namespace io
