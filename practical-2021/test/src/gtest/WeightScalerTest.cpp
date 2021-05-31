#include "gtest/gtest.h"
#include "../../../src/WeightScaler.hpp"

#include <random>
#include <iomanip>


class WeightScalerTest : public testing::Test {
protected:
  void checkDuplicates(const std::vector<std::vector<int64_t>>& scaled_weights, const std::vector<std::vector<double>>& weights){
    for(size_t i_1 = 0; i_1 < weights.size(); ++i_1){
      for(size_t j_1 = 0; j_1 < weights.size(); ++j_1){
        for(size_t i_2 = 0; i_2 < weights.size(); ++i_2){
          for(size_t j_2 = 0; j_2 < weights.size(); ++j_2){
            if (i_1 == i_2 && j_1 == j_2) continue;
            EXPECT_TRUE(scaled_weights[i_1][j_1] != scaled_weights[i_2][j_2]);
            if (scaled_weights[i_1][j_1] == scaled_weights[i_2][j_2]) {
              //std::cout << weights[i_1][j_1] << " " << weights[i_2][j_2] << std::endl;
              //std::cout << scaled_weights[i_1][j_1] << " " << scaled_weights[i_2][j_2] << std::endl;
            }
          }
        }
      }
    }
  }
};


//TODO move randomness out of here
TEST_F(WeightScalerTest, test_maximum) {
  size_t n = 10;
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distr(-10, 10);
  std::vector<std::vector<double>> weights = std::vector<std::vector<double>> (n, std::vector<double>(n));
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j < n; ++j){
      weights[i][j] = distr(eng);
    }
  }
  std::vector<std::vector<int64_t>> scaled_weights = WeightScaler::scale(weights);
  checkDuplicates(scaled_weights, weights);


}
