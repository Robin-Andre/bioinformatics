#pragma once
#include <iostream>
#include <vector>
#include "ortools/graph/assignment.h"
#include <utility>
#include "WeightScaler.hpp"

class MaximumMatcher {
public:
  static int convert_weight(double weight) {
    return (int) (-256 * (weight + 0.5));
  }
  static std::vector<size_t> match_vector(const std::vector<std::vector<double>>& weights) {
    std::vector<size_t> matching_vector(weights.size());
    operations_research::SimpleLinearSumAssignment assignment;
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        //TODO the transformation between probabilities and whole numbers recheck. Also is this the integer rounding?
        //assignment.AddArcWithCost(i, j, (max_weight.first - weights[i][j])/min_diff);
        assignment.AddArcWithCost(i, j,  convert_weight(weights[i][j]));
      }
    }
    double result = 0;
    if (assignment.Solve() == operations_research::SimpleLinearSumAssignment::OPTIMAL) {
      //printf("A perfect matching exists.\n");
      //printf("The best possible cost is %d.\n", assignment.OptimalCost());
      //printf("An optimal assignment is:\n");
      for (int node = 0; node < assignment.NumNodes(); ++node) {

        /*printf("left node %d assigned to right node %d with cost %i (%i).\n",
        node,
        assignment.RightMate(node),
        (- assignment.AssignmentCost(node)),
        assignment.AssignmentCost(node));*/
        matching_vector[node] = assignment.RightMate(node);
        //std::cout << "Result(Node): " << node << "->" << assignment.RightMate(node) << " = " << result << "\n";
      }
      //printf("Note that it may not be the unique optimal assignment.");
    }
    return matching_vector;
  }
  static double match(const std::vector<std::vector<double>>& weights) {
    std::vector<size_t> matching = match_vector(weights);
    double result = 0;
    for(size_t i = 0; i < matching.size(); ++i) {
      result += weights[i][matching[i]];
      std::cout << "Result: " << result << "\n";
    }
    return result;
  }
private:
  static std::pair<double,double> findMax(const std::vector<std::vector<double>>& weights) {
    std::pair<double, double> max_weights = std::make_pair(-DBL_MAX, -DBL_MAX);
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        if(weights[i][j] > max_weights.first){
          max_weights.second = max_weights.first;
          max_weights.first = weights[i][j];
        }
      }
    }
    return max_weights;
  }

};
