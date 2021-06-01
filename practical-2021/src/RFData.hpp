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

class RFData {
public:
  size_t tree_count;
  size_t unique_count;
  size_t tip_count;
  std::vector<double> distances;
  double average_distance;

  RFData(){};
  RFData(size_t tree_count, size_t tip_count) : tree_count(tree_count), tip_count(tip_count){};

  std::vector<double> getRelativeDistances() const {
    std::vector<double> relative_distances;
    for(size_t i = 0; i < distances.size(); ++i){
      relative_distances.emplace_back((double) distances[i] / (2*(tip_count-3)));
    }
    return relative_distances;
  }

  //tip count is not compared, as not always given... maybe needs to be added
  bool operator==(const RFData& other) const {
    bool is_eq = distancesEqual(other.distances);
    is_eq &= near(average_distance, other.average_distance);
    is_eq &= (unique_count == other.unique_count);
    is_eq &= (tree_count == other.tree_count);
    return is_eq;
  }

private:
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


};
