#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../datastructures/PllTree.hpp"
#include "CommandLineOptions.hpp"
#include "IOData.h"
#include "../RFData.hpp"

enum Format {RAXML};
class RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data);

};

class RAXMLWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data) {
    std::ofstream out_file(path + "/distances");
    size_t tree_count = data.getTreeCount();
    if (out_file.is_open()) {
      std::vector<double> relative_distances = data.getRelativeDistances();
      std::vector<double> distances = data.getDistances();
      size_t k = 0;
      for(unsigned i = 0; i < tree_count; ++i) {
        for(unsigned j = i + 1; j < tree_count; ++j) {
          out_file << i << " " << j << " " << distances[k] << " " << relative_distances[k] << std::endl;
          ++k;
        }
      }
      out_file.close();
    } else {
      throw ("Cannot write to " + path + "/distances");
    }
    out_file = std::ofstream(path + "/info");
    if (out_file.is_open()) {
      out_file << "Found " << tree_count << " trees in File" << std::endl;
      out_file << "Number of unique trees in this tree set: " << data.getUniqueCount() << std::endl;
      out_file << "Average relative RF in this set: " << data.getAverageDistance() << std::endl;
      out_file.close();
    } else {
      throw ("Cannot write to " + path +  "/info");
    }
  }

};


class JSONWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data) {
    std::ofstream out_stream(path + "/results.json");
    if (out_stream.is_open()) {
        nlohmann::json j;
        io::to_json(j, data.getIOData());
        j.dump(4);
        out_stream << j;
        out_stream.close();
    } else {
      throw ("Cannot write JSON to " + path + "/results.json");
    }
  }

};


class MatrixWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data) {
    size_t tree_count = data.getTreeCount();
    std::ofstream out_file(path);
    if (out_file.is_open()) {
      std::vector<std::vector<double>> matrix = data.getFullMatrix();
      unsigned k = 0;
      for(unsigned i = 0; i < tree_count; ++i) {
        for(unsigned j = 0; j < tree_count - 1; ++j) {
          out_file << matrix[i][j]  << " ";
          ++k;
        }
        out_file << matrix[i][tree_count-1] << std::endl;
      }
      out_file.close();
    } else {
      throw ("Cannot write to " + path);
    }
  }


};
