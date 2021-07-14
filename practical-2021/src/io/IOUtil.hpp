#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include <vector>
#include <iostream>
#include <numeric>

class IOUtil {
public:


static std::vector<std::vector<double>> fullMatrixToHalfMatrix(const std::vector<std::vector<double>>& full) {
  std::vector<std::vector<double>> half(full.size(), std::vector<double>());
  for(size_t i = 0; i < full.size(); ++i){
    assert(full[i].size() == full.size());
    for(size_t j = i; j < full.size(); ++j){
      //seems as if there is not more precision in input :o
      //assert(std::abs(full[i][j] - full[j][i]) < 0.001);
      half[j].emplace_back(full[j][i]);
    }
  }
  return half;
}


static std::vector<std::vector<double>> halfMatrixToFullMatrix(const std::vector<std::vector<double>>& half) {
  std::vector<std::vector<double>> full;
  for(size_t i = 0; i < half.size(); ++i){
    std::vector<double> row(half.size());
    for(size_t j = 0; j < i+1; ++j){
      row[j] = half[i][j];
    }
    for(size_t j = 0; j < half.size()-i-1; ++j){
      row[i+j+1] = half[i+j+1][i];
    }
    full.emplace_back(row);
  }
  return full;
}

static size_t calculateUniqueCount(const std::vector<std::vector<double>>& matrix) {
  size_t unique_count = matrix.size();
  for(size_t i = 0; i < matrix.size(); ++i){
    bool is_unique = true;
    for(size_t j = 0; j < i; ++j){
      if (matrix[i][j] <= 0 && is_unique){
        is_unique = false;
        --unique_count;
      }
    }
  }
  assert(unique_count >= 1);
  return unique_count;
}

static std::vector<std::vector<double>> vectorToHalfMatrix(const std::vector<double>& distances, size_t tree_count) {
  std::vector<std::vector<double>> half(tree_count, std::vector<double>());
  size_t k = 0;
  for(size_t i = 0; i < tree_count; ++i){
    for(size_t j = i+1; j < tree_count; ++j){
      half[j].emplace_back(distances[k]);
      ++k;
    }
  }
  //fill diagonal
  for(size_t i = 0; i < tree_count; ++i){
    half[i].emplace_back(0);
  }
  return half;
}

static double calculateAverageDistance(const std::vector<std::vector<double>>& half) {
  double result = 0;
  size_t k = 0;
  for(size_t i = 0; i < half.size(); ++i) {
    for(size_t j = i+1; j < half.size(); ++j) {
      result += half[j][i];
      ++k;
    }
  }
  return result / static_cast<double>(k);
}


};
