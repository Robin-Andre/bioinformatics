#include "gtest/gtest.h"
#include "io/IOUtil.hpp"
#include <stdlib.h>
#include <random>
#include <iomanip>
#include <vector>
class IOUtilTest : public testing::Test {

};


/**
* Test functions of IOUtils using data structures filled with random values
*
*/

TEST_F(IOUtilTest, test_half_to_full){
  size_t n = 10;
  std::vector<std::vector<double>> half = std::vector<std::vector<double>>(n, std::vector<double>());
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distr(-10, 10);
  for(size_t i = 0; i < n; ++i){
    for(size_t j = i; j < n; ++j){
      half[j].emplace_back(distr(eng));
    }
  }
  std::vector<std::vector<double>> full = IOUtil::halfMatrixToFullMatrix(half);
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j <= i; ++j){
      EXPECT_EQ(full[i][j], half[i][j]);
      EXPECT_EQ(full[i][j], full[j][i]);
    }
  }
}

TEST_F(IOUtilTest, test_full_to_half){
  size_t n = 10;
  std::vector<std::vector<double>> full = std::vector<std::vector<double>>(n, std::vector<double>(n));
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distr(-10, 10);
  for(size_t i = 0; i < n; ++i){
    full[i][i] = distr(eng);
    for(size_t j = i+1; j < n; ++j){
      double r = distr(eng);
      full[i][j] = r;
      full[j][i] = r;
    }
  }
  std::vector<std::vector<double>> half = IOUtil::fullMatrixToHalfMatrix(full);
  for(size_t i = 0; i < n; ++i){
    for(size_t j = 0; j <= i; ++j){
      EXPECT_EQ(full[i][j], half[i][j]);
    }
  }
}


TEST_F(IOUtilTest, test_vector_to_half){
  size_t n = 10;
  std::vector<double> vec;
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> distr(-10, 10);
  for(size_t i = 0; i < n; ++i){
    for(size_t j = i+1; j < n; ++j){
      vec.emplace_back(distr(eng));
    }
  }
  std::vector<std::vector<double>> half = IOUtil::vectorToHalfMatrix(vec, n);
  size_t k = 0;
  for(size_t i = 0; i < n; ++i) {
    for(size_t j = i+1; j < n; ++j) {
      EXPECT_EQ(half[j][i], vec[k]);
      ++k;
    }
  }
}
