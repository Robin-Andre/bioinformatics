#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../PllTree.hpp"
#include "../RFDistance.hpp"
#include "CommandLineOptions.hpp"

enum Format {RAXML};
class RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data);

};

class RAXMLWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data) {
    std::ofstream out_file(path);
    unsigned k = 0;
    for(unsigned i = 0; i < data.tree_count; ++i) {
      for(unsigned j = i + 1; j < data.tree_count; ++j) {
        out_file << i << " " << j << " " << data.distances[k] << " " << data.relative_distances[k] << std::endl;
        ++k;
      }
    }
    out_file.close();
    std::cout << "Result File written at: "<< path << std::endl;
    std::cout << "Unique count: " << data.unique_count
              << " Average distance: " << data.average_distance
              << " Tree count: " << data.tree_count << "\n";
  }

private:
};
