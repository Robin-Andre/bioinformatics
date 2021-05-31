#pragma once
#include <iostream>
#include <vector>
#include "ortools/graph/assignment.h"
#include <utility>

class MaximumMatcher {
public:
  static double match(const std::vector<std::vector<double>>& weights) {
    operations_research::SimpleLinearSumAssignment assignment;
    std::pair<double,double> max_weight = findMax(weights);
    double min_diff = max_weight.first - max_weight.second; //This has an issue
    int larger_than_max_weight_int = (int) (2 * std::max(max_weight.first, max_weight.second));
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        //TODO the transformation between probabilities and whole numbers recheck. Also is this the integer rounding?
        //assignment.AddArcWithCost(i, j, (max_weight.first - weights[i][j])/min_diff);
        assignment.AddArcWithCost(i, j,  -(4 * weights[i][j] + 0.5));
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
        result += weights[node][assignment.RightMate(node)];
        //std::cout << "Result(Node): " << node << "->" << assignment.RightMate(node) << " = " << result << "\n";
      }
      //printf("Note that it may not be the unique optimal assignment.");
    } else {
      //printf("There is an issue with the input or no perfect matching exists.");
      return -DBL_MAX;
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
