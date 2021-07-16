#pragma once
#include <iostream>
#include <vector>
#include "ortools/graph/assignment.h"
#include <utility>

/**
 * Offers functionality to determine a maximum matching for weights
 * corresponding to edge weights in a bipartite graph
 */
class MaximumMatcher {
public:
  /**
   * Determines a maximum matching for a bipartite graph with the provided weights
   *
   * @param weights: edge weights for the bipartite graph
   * @return permutation coresponding a maximum matching
   */

  static std::vector<size_t> match_vector(const std::vector<std::vector<double>>& weights) {
    std::vector<size_t> matching_vector(weights.size());
    operations_research::SimpleLinearSumAssignment assignment;
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        assignment.AddArcWithCost(static_cast<int>(i), static_cast<int>(j),
          convert_weight(weights[i][j])); //@Softwipe, added explicit conversion
      }
    }
    if (assignment.Solve() == operations_research::SimpleLinearSumAssignment::OPTIMAL) {
      size_t num_nodes = static_cast<size_t> (assignment.NumNodes());
      for (size_t node = 0; node < num_nodes; ++node) {
        matching_vector[node] = static_cast<size_t>(assignment.RightMate(static_cast<int>(node)));
        //std::cout << "Result(Node): " << node << "->" << assignment.RightMate(node) << " = " << result << "\n";
      }
    }
    return matching_vector;
  }

  /**
   * Determines the weight a maximum matching for a bipartite
   * graph with the provided weights
   *
   * @param weights: edge weights for the bipartite graph
   * @return weight of a maximum matching
   */

  static double match(const std::vector<std::vector<double>>& weights) {
    std::vector<size_t> matching = match_vector(weights);
    assert(matching.size() == weights.size());
    double result = 0;
    for(size_t i = 0; i < matching.size(); ++i) {
      result += weights[i][matching[i]];
    }
    return result;
  }

private:
  /**
   * Rescales the weights to make them proper for ortools linear solver
   *
   * @param weights: weights to be converted
   * @return converted weights
   */

  static int convert_weight(double weight) {
    return static_cast<int> ((-65536 * weight) + 0.5);
  }

};
