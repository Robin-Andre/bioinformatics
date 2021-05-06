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
    std::cout << "file created" << std::endl;
    std::ofstream out_file(path + "/distances");
    std::cout << "file created" << std::endl;
    unsigned k = 0;
    for(unsigned i = 0; i < data.tree_count; ++i) {
      for(unsigned j = i + 1; j < data.tree_count; ++j) {
        out_file << i << " " << j << " " << data.distances[k] << " " << data.relative_distances[k] << std::endl;
        ++k;
      }
    }

    out_file = std::ofstream(path + "/info");
    out_file << "Found " << data.tree_count << " trees in File" << std::endl;
    out_file << "Number of unique trees in this tree set: " << data.unique_count << std::endl;
    out_file << "Average relative RF in this set: " << data.average_distance << std::endl;
    out_file.close();
  }

};
