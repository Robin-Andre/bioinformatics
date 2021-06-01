#pragma once
#include<vector>
#include <fstream>
#include <iostream>
#include "../datastructures/PllTree.hpp"
#include "../RFDistance.hpp"
#include "CommandLineOptions.hpp"
#include "IOData.h"

enum Format {RAXML};
class RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data);

};

class RAXMLWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data) {
    std::ofstream out_file(path + "/distances");
    if (out_file.is_open()) {
      std::vector<double> relative_distances = data.getRelativeDistances();
      unsigned k = 0;
      for(unsigned i = 0; i < data.tree_count; ++i) {
        for(unsigned j = i + 1; j < data.tree_count; ++j) {
          out_file << i << " " << j << " " << data.distances[k] << " " << relative_distances[k] << std::endl;
          ++k;
        }
      }
      out_file.close();
    } else {
      throw ("Cannot write to " + path + "/distances");
    }
    out_file = std::ofstream(path + "/info");
    if (out_file.is_open()) {
      out_file << "Found " << data.tree_count << " trees in File" << std::endl;
      out_file << "Number of unique trees in this tree set: " << data.unique_count << std::endl;
      out_file << "Average relative RF in this set: " << data.average_distance << std::endl;
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
        io::to_json(j, convertToIOData(data));
        j.dump(4);
        out_stream << j;
        out_stream.close();
    } else {
      throw ("Cannot write JSON to " + path + "/results.json");
    }
  }

private:

  static io::IOData convertToIOData(const RFData& rf_data){
    io::IOData io_data;
    io_data.mean_rf_dst = rf_data.average_distance;
    io_data.split_score_calc = io::IOData::Metric::RF;
    io_data.mean_modified_rf_dst = 0;
    for(size_t i = 0; i < rf_data.tree_count; i++){
      io_data.taxa_names.emplace_back("taxa_" + std::to_string(i));
    }
    io_data.pairwise_distance_mtx = toMatrix(rf_data.distances, rf_data.tree_count);
    io_data.git_revision = "insert git revision here";
    io_data.cpuInformation = "insert cpu information here";
    io_data.number_of_unique_trees = rf_data.unique_count;
    return io_data;
  }

  static std::vector<std::vector<double>> toMatrix(const std::vector<double>& distances, size_t tree_count){
    std::vector<std::vector<double>> matrix;
    size_t k = 0;
    for(size_t i = 0; i < tree_count-1; ++i){
      std::vector<double> row;
      for(size_t j = i+1; j < tree_count; ++j){
        row.emplace_back(distances[k]);
        ++k;
      }
      matrix.emplace_back(row);
    }
    return matrix;
  }

};


class MatrixWriter : public RFDataWriter {
public:
  static void write(const std::string& path, const RFData& data) {
    std::ofstream out_file(path + "/distances");
    if (out_file.is_open()) {
      std::vector<std::vector<double>> matrix = vectorToMatrix(data.distances, data.tree_count);
      unsigned k = 0;
      for(unsigned i = 0; i < data.tree_count; ++i) {
        for(unsigned j = 0; j < data.tree_count - 1; ++j) {
          out_file << matrix[i][j]  << " ";
          ++k;
        }
        out_file << matrix[i][data.tree_count-1] << std::endl;
      }
      out_file.close();
    } else {
      throw ("Cannot write to " + path + "/distances");
    }
  }
private:
  static std::vector<std::vector<double>> vectorToMatrix(std::vector<double> v, size_t n) {
    std::vector<std::vector<double>> matrix = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
    size_t k  = 0;
    for(size_t i = 0; i < n; ++i){
      for(size_t j = i+1; j < n; ++j){
        matrix[i][j] = v[k];
        matrix[j][i] = v[k];
        ++k;
      }
    }
    return matrix;
  }

};
