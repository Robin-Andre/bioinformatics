#pragma once
#include <iostream>
#include <vector>
#include "ortools/graph/assignment.h"

class MaximumMatcher {
public:
  static double match(const std::vector<std::vector<double>>& weights) {
    operations_research::SimpleLinearSumAssignment assignment;
    double max_weight = findMax(weights);
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        assignment.AddArcWithCost(i, j, max_weight - weights[i][j]);
      }
    }
    double result = 0;
    if (assignment.Solve() == operations_research::SimpleLinearSumAssignment::OPTIMAL) {
      //printf("A perfect matching exists.\n");
      //printf("The best possible cost is %d.\n", assignment.OptimalCost());
      //printf("An optimal assignment is:\n");
      for (int node = 0; node < assignment.NumNodes(); ++node) {
        /*printf("left node %d assigned to right node %d with cost %f (%f).\n",
        node,
        assignment.RightMate(node),
        (max_weight - assignment.AssignmentCost(node)),
        assignment.AssignmentCost(node));*/
        result += (max_weight - assignment.AssignmentCost(node));
      }
      //printf("Note that it may not be the unique optimal assignment.");
    } else {
      //printf("There is an issue with the input or no perfect matching exists.");
      return -DBL_MAX;
    }
    return result;
  }
private:
  static double findMax(const std::vector<std::vector<double>>& weights) {
    double max_weight = -DBL_MAX;
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        max_weight = std::max(max_weight, weights[i][j]);
      }
    }
    return max_weight;
  }

};
