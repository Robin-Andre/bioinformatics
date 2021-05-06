#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../PllTree.hpp"
#include "../RFDistance.hpp"
#include "CommandLineOptions.hpp"
namespace io {
/*static std::vector<PllTree> readTreeFile(const std::string& filepath) {
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
}*/

enum Format {RAXML};

static void writeRAXML(const RFData& result, const std::string& path) {
  std::ofstream out_file(path);
  unsigned k = 0;
  for(unsigned i = 0; i < result.tree_count; ++i) {
    for(unsigned j = i + 1; j < result.tree_count; ++j) {
      out_file << i << " " << j << " " << result.distances[k] << " " << result.relative_distances[k] << std::endl;
      ++k;
    }
  }
  out_file.close();
}

static void writeRFData(const RFData& result, const Config& config) {
  //Maybe an argument later
  Format format = RAXML;

  if(config.output_file_path.size() <= 1) {
    std::cerr << "The Output was not properly specified...exiting";
    exit(1);
  }
  switch (format) {
    case RAXML:
    {
      writeRAXML(result, config.output_file_path);
      break; //@Luise <- This one was missing :P
    }
    default:
      throw ("No proper output format specified!");
  }
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



static size_t getAbsoluteDistanceFromString(const std::string &line) {
    std::istringstream iss (line);
    std::string item;
    size_t i = 0;
    while (std::getline(iss, item, ' ') && i < 2) {
        i++;
    }
    return std::stoi(item);
}

static float getRelativeDistanceFromString(const std::string &line) {
    std::istringstream iss (line);
    std::string item;
    size_t i = 0;
    while (std::getline(iss, item, ' ') && i < 3) {
        i++;
    }
    return std::stof(item);
}


static std::string parseAfterPrefix(const std::string& line, const std::string& prefix){
  std::string result;
  std::string dummy;
  std::istringstream iss (line);
  std::istringstream iss2 (prefix);
  while (std::getline(iss2, dummy, ' ')){
    std::getline(iss, result, ' ');
  }
  std::getline(iss, result, ' ');
  return result;

}



static RFData readRAXML(const std::string& path) {
  RFData data;
  std::fstream res_file;
  res_file.open(path + "/RAxML_RF-Distances.0"  ,std::ios::in);
  if (res_file.is_open()){
    std::string line;
    std::vector<std::string> parts;
    while(std::getline(res_file, line)){
      data.distances.push_back(getAbsoluteDistanceFromString(line));
      data.relative_distances.push_back(getRelativeDistanceFromString(line));
    }
    res_file.close(); //close the file object.
  } else {
    throw (path +  "/RAxML_RF-Distances.0 not found!");
  }


  std::vector<std::string> prefixes{"Found ",
   "Number of unique trees in this tree set: ",
   "Average relative RF in this set: "};

  res_file.open(path + "/RAxML_info.0"  ,std::ios::in);
  std::vector<std::string> prefix_matches(prefixes.size());
  if (res_file.is_open()){
    std::string line;
    while(std::getline(res_file, line)){
      for(size_t i = 0; i < prefixes.size(); ++i){
        auto prefix = prefixes[i];
        auto match = std::mismatch(prefix.begin(), prefix.end(), line.begin());
        if (match.first == prefix.end())
        {
          prefix_matches[i] = parseAfterPrefix(line, prefix);
        }
      }

    }
    res_file.close();
  } else {
    throw (path +  "/RAxML_info.0 not found!");
  }
  data.tree_count = std::stoi(prefix_matches[0]);
  data.unique_count = std::stoi(prefix_matches[1]);
  data.average_distance = std::stof(prefix_matches[2]);
  return data;
}

static RFData readRFData(const std::string& path, Format format) {
  switch(format) {
    case RAXML:
    {
      return readRAXML(path);
    }
    default:
      throw "No proper input format specified!";
    }
}


}//namespace io
