#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllTree.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <cfloat>
#include <numeric>
#include "io/IOData.h"

class RFData {
public:

  RFData(size_t tree_count, size_t tip_count, size_t unique_count, const std::vector<double>& distances, bool relative_average) :
  tree_count(tree_count),
  tip_count(tip_count),
  unique_count(unique_count),
  distances(distances){
    assert(distances.size() == (tree_count * (tree_count - 1)) / 2);
    average_distance = calculateAverageDistance(relative_average);
  };

  RFData(size_t tree_count, size_t unique_count, double average_distance, const std::vector<double>& distances) :
  tree_count(tree_count),
  tip_count(-1),
  unique_count(unique_count),
  average_distance(average_distance),
  distances(distances){};

  RFData(const std::vector<std::vector<double>>& matrix){
    tree_count = matrix.size();
    tip_count = -1;
    unique_count = calculateUniqueCount(matrix);
    distances = fullMatrixToVector(matrix);
    average_distance = calculateAverageDistance(false);
  }

  RFData(const io::IOData& io_data){
    average_distance = io_data.mean_rf_dst;
    tree_count = io_data.pairwise_distance_mtx.size() + 1;
    unique_count = io_data.number_of_unique_trees;
    distances = halfMatrixToVector(io_data.pairwise_distance_mtx);
  }

  std::vector<double> getRelativeDistances() const {
    std::vector<double> relative_distances;
    for(size_t i = 0; i < distances.size(); ++i){
      relative_distances.emplace_back((double) distances[i] / (2*(tip_count-3)));
    }
    return relative_distances;
  }

 std::vector<std::vector<double>> getFullMatrix() const{
    std::vector<std::vector<double>> matrix = std::vector<std::vector<double>>(tree_count, std::vector<double>(tree_count, 0));
    size_t k  = 0;
    for(size_t i = 0; i < tree_count; ++i){
      for(size_t j = i+1; j < tree_count; ++j){
        matrix[i][j] = distances[k];
        matrix[j][i] = distances[k];
        ++k;
      }
    }
    return matrix;
  }

  io::IOData getIOData()const{
    io::IOData io_data;
    io_data.mean_rf_dst = average_distance;
    io_data.split_score_calc = io::IOData::Metric::RF;
    io_data.mean_modified_rf_dst = 0;
    for(size_t i = 0; i < tree_count; i++){
      io_data.taxa_names.emplace_back("taxa_" + std::to_string(i));
    }
    io_data.pairwise_distance_mtx = getHalfMatrix();
    io_data.git_revision = "insert git revision here";
    io_data.cpuInformation = "insert cpu information here";
    io_data.number_of_unique_trees = unique_count;
    return io_data;
  }



  //tip count is not compared, as not always given... maybe needs to be added
  bool operator==(const RFData& other) const {
    bool is_eq = distancesEqual(other.distances);
    is_eq &= near(average_distance, other.average_distance);
    is_eq &= (unique_count == other.unique_count);
    is_eq &= (tree_count == other.tree_count);
    return is_eq;
  }

  std::vector<double> getDistances() const {return distances;}
  size_t getTreeCount() const {return tree_count;}
  size_t getUniqueCount() const {return unique_count;}
  double getAverageDistance() const {return average_distance;}

private:

  size_t tree_count;
  size_t unique_count;
  size_t tip_count;
  std::vector<double> distances;
  double average_distance;

bool near(double a, double b) const {
  auto absA = std::abs(a);
  auto absB = std::abs(b);
  auto largest = (absA < absB) ? absB : absA;
  auto smallest = (absA < absB) ? absA : absB;
  return largest - smallest <= largest * static_cast<double>(FLT_EPSILON) * 1e5;
}

bool distancesEqual(std::vector<double> otherDistances) const{
  if(distances.size() != otherDistances.size()) return false;
  for(size_t i = 0; i < distances.size(); ++i){
    if(! near(distances[i], otherDistances[i])) return false;
  }
  return true;
}


std::vector<double> fullMatrixToVector(const std::vector<std::vector<double>>& matrix) const{
  std::vector<double> result;
  for(size_t i = 0; i < matrix.size(); ++i){
    assert(matrix[i].size() ==  matrix.size());
    for(size_t j = i+1; j < matrix.size(); ++j){
      result.emplace_back(matrix[i][j]);
    }
  }
  return result;
}

std::vector<double> halfMatrixToVector(const std::vector<std::vector<double>>& matrix) {
  std::vector<double> result;
  for(std::vector<double> row : matrix)  {
    for(double el : row){
      result.emplace_back(el);
    }
  }
  return result;
}


size_t calculateUniqueCount(const std::vector<std::vector<double>>& matrix) const{
  size_t unique_count = matrix.size();
  bool is_unique = true;
  for(size_t i = 0; i < matrix.size(); ++i){
    for(size_t j = i+1; j < matrix.size(); ++j){
      if (matrix[i][j] == 0 && is_unique){
        is_unique = false;
        --unique_count;
      }
    }
  }
  return unique_count;
}

std::vector<std::vector<double>> getHalfMatrix() const{
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

double calculateAverageDistance(bool relative) {
  if (relative) {
    return calculateAverageDistance(false) / (2*(tip_count-3));
  } else {
    return (std::accumulate(distances.begin(), distances.end(), 0.0)) / distances.size();
  }
}


};
