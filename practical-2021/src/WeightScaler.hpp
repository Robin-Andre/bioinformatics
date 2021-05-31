#pragma once
#include <vector>

class WeightScaler {
public:
  static std::vector<std::vector<int64_t>> scale(const std::vector<std::vector<double>>& weights) {
    double high_value = 100;
    std::vector<std::vector<int64_t>>  result = std::vector<std::vector<int64_t>> (weights.size(), std::vector<int64_t>(weights.size()));
    double min_diff = findMinDiff(weights);
    for (size_t i = 0; i < weights.size(); ++i) {
      for(size_t j = 0; j < weights[i].size(); ++j) {
        result[i][j] = (int64_t)((high_value - weights[i][j])/min_diff);
        //std::cout << weights[i][j] << " >> " << result[i][j] << std::endl;
      }
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


  static double findMinDiff(const std::vector<std::vector<double>>& weights){
    double min_diff = DBL_MAX;
    for(size_t i_1 = 0; i_1 < weights.size(); ++i_1){
      for(size_t j_1 = 0; j_1 < weights.size(); ++j_1){
        for(size_t i_2 = 0; i_2 < weights.size(); ++i_2){
          for(size_t j_2 = 0; j_2 < weights.size(); ++j_2){
            if (i_1 == i_2 && j_1 == j_2) continue;
            double diff = std::abs(weights[i_1][j_1]-  weights[i_2][j_2]);
            if (diff != 0) min_diff = std::min(min_diff, diff);
            if(weights[i_1][j_1] != 0) min_diff = std::min(min_diff, std::abs(weights[i_1][j_1]));
          }
        }
      }
    }
  }

};
