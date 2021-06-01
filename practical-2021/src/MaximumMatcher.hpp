#pragma once
#include <iostream>
#include <vector>
#include "ortools/graph/assignment.h"
#include <utility>

class MaximumMatcher {
public:
  static double match(const std::vector<std::vector<double>>& weights) {
    operations_research::SimpleLinearSumAssignment assignment;
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        assignment.AddArcWithCost(i, j, scale(weights[i][j]));
      }
    }
    double result = 0;
    if (assignment.Solve() == operations_research::SimpleLinearSumAssignment::OPTIMAL) {
      //printf("A perfect matching exists.\n");
      //printf("The best possible cost is %d.\n", assignment.OptimalCost());
      //printf("An optimal assignment is:\n");
      for (int node = 0; node < assignment.NumNodes(); ++node) {
        /*printf("left node %d assigned to right node %d with cost %f.\n",
        node,
        assignment.RightMate(node),
        assignment.AssignmentCost(node));*/
        result += weights[node][assignment.RightMate(node)];
      }
      //printf("Note that it may not be the unique optimal assignment.");
    } else {
      //printf("There is an issue with the input or no perfect matching exists.");
      return -DBL_MAX;
    }
    return result;
  }
private:
  static size_t scale(double weight){
    double multiplicator = 10000000;
    return std::round(-multiplicator*weight);
  }

};
