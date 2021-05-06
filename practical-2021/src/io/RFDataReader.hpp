#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../PllTree.hpp"
#include "../RFDistance.hpp"
#include "CommandLineOptions.hpp"

class RFDataReader {
public:
  static RFData read(const std::string& path);

};

class RAXMLReader : public RFDataReader {
public:
  static RFData read(const std::string& path) {
    RFData data;
    std::fstream res_file;
    res_file.open(path + "/distances"  ,std::ios::in);
    if (res_file.is_open()){
      std::string line;
      std::vector<std::string> parts;
      while(std::getline(res_file, line)){
        data.distances.push_back(getAbsoluteDistanceFromString(line));
        data.relative_distances.push_back(getRelativeDistanceFromString(line));
      }
      res_file.close(); //close the file object.
    } else {
      throw (path +  "/distances not found!");
    }


    std::vector<std::string> prefixes{"Found ",
     "Number of unique trees in this tree set: ",
     "Average relative RF in this set: "};

    res_file.open(path + "/info"  ,std::ios::in);
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
      throw (path +  "/info not found!");
    }
    data.tree_count = std::stoi(prefix_matches[0]);
    data.unique_count = std::stoi(prefix_matches[1]);
    data.average_distance = std::stof(prefix_matches[2]);
    return data;
  }

private:
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

};
