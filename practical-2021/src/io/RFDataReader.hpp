#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../datastructures/PllTree.hpp"
#include "CommandLineOptions.hpp"
#include "IOData.h"
#include "../RFData.hpp"

class RFDataReader {
public:
  static RFData read(const std::string& path);

};

class RAXMLReader : public RFDataReader {
public:
  static RFData read(const std::string& path) {
    std::vector<double> distances;
    std::fstream res_file;
    res_file.open(path + "/distances"  ,std::ios::in);
    if (res_file.is_open()){
      std::string line;
      while(std::getline(res_file, line)){
        distances.push_back(getAbsoluteDistanceFromString(line));
        //data.relative_distances.push_back(getRelativeDistanceFromString(line));
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
    return RFData(std::stoi(prefix_matches[0]), std::stoi(prefix_matches[1]), std::stod(prefix_matches[2]), distances);
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


class JSONReader : public RFDataReader {
public:
  static RFData read(const std::string& path) {
    io::IOData io_data;
    std::fstream res_file;
    res_file.open(path + "/results.json", std::ios::in);
    if (res_file.is_open()){
      std::string json_str;
      std::string line;
      while(std::getline(res_file, line)){
        json_str += line;
        json_str += '\n';
      }
      nlohmann::json jsonIn = nlohmann::json::parse(json_str);
      io::from_json(jsonIn, io_data);
      return RFData(io_data);
    } else {
      throw ("Cannot read JSON from " + path  + "/results.json");
    }

  }

};


class MatrixReader : RFDataReader {
public:
  static RFData read(const std::string& path) {
    std::fstream res_file;
    res_file.open(path, std::ios::in);
    if (res_file.is_open()){
      std::vector<std::vector<double>> matrix;
      std::string line;
      while(std::getline(res_file, line)){
        std::vector<double> row;
        std::string token;
        size_t pos = 0;
        while ((pos = line.find(" ")) != std::string::npos) {
            row.push_back(std::stold(line.substr(0, pos)));
            line.erase(0, pos + 1);
        }
        row.push_back(std::stof(line));
        matrix.push_back(row);
      }
      res_file.close();
      return RFData(matrix);
    } else {
      throw (path +  " not found!");
    }

  }

};
