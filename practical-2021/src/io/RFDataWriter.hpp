#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../datastructures/PllTree.hpp"
#include "CommandLineOptions.hpp"
#include "IOData.hpp"
#include "IOUtil.hpp"

enum Format {RAXML};

/*
 * Abstract class for writing IOData to file
 */
class RFDataWriter {
public:
  static void write(const std::string& path, const io::IOData& data);

};

/*
 * Class for writing IOData to file in the format used by RAXML
 */
class RAXMLWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const io::IOData& data) {
    std::ofstream out_file(path + "/distances");
    size_t tree_count = data.pairwise_distance_mtx.size();
    if (out_file.is_open()) {
      for(size_t i = 0; i < tree_count; ++i) {
        for(size_t j = i+1; j < tree_count; ++j) {
          assert(data.pairwise_distance_mtx[j].size() > i);
          out_file << i << " " << j << " " << data.pairwise_distance_mtx[j][i] << std::endl;
        }
      }
      out_file.close();
    } else {
      std::cerr << "Cannot write to " << path << "/distances\n";
      exit(1);
    }
    out_file = std::ofstream(path + "/info");
    if (out_file.is_open()) {
      out_file << "Found " << tree_count << " trees in File" << std::endl;
      out_file << "Number of unique trees in this tree set: " << data.number_of_unique_trees << std::endl;
      out_file << "Average relative RF in this set: " << data.mean_dst << std::endl;
      out_file.close();
    } else {
      std::cerr << "Cannot write to " << path << "/info\n";
      exit(1);
    }
  }

};

/*
 * Class for writing IOData to file in JSON Format
 */
class JSONWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const io::IOData& data) {
    std::ofstream out_stream(path + ".json");
    if (out_stream.is_open()) {
        nlohmann::json j;
        io::to_json(j, data);
        j.dump(4);
        out_stream << j;
        out_stream.close();
    } else {
      std::cerr << "Cannot write JSON to " << path << "/.json\n";
      exit(1);
    }
  }

};


/*
 * Class for writing IOData to file in Matrix Format
 */
class MatrixWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const io::IOData& data) {
    size_t tree_count = data.pairwise_distance_mtx.size();
    std::ofstream out_file(path);
    if (out_file.is_open()) {
      std::vector<std::vector<double>> matrix = IOUtil::halfMatrixToFullMatrix(data.pairwise_distance_mtx);
      for(size_t i = 0; i < tree_count; ++i) {
        assert(matrix[i].size() == tree_count);
        for(size_t j = 0; j < tree_count - 1; ++j) {
          out_file << matrix[i][j]  << " ";
        }
        out_file << matrix[i][tree_count-1] << std::endl;
      }
      out_file.close();
    } else {
      std::cerr << "Cannot write to " << path << "\n";
      exit(1);
    }
  }


};
