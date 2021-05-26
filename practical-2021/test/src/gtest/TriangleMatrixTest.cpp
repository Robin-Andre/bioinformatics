#include "gtest/gtest.h"
#include "../../../src/datastructures/TriangleMatrix.hpp"
#include <stdlib.h>
class TraingleMatrixTest : public testing::Test {

};

TEST_F(TraingleMatrixTest, simple_test_diagonal){
  size_t n = 10;
  std::vector<std::vector<size_t>> reference_matrix = std::vector<std::vector<size_t>>(n, std::vector<size_t>(n));
  TriangleMatrix<size_t> triangle_matrix = TriangleMatrix<size_t>(n, true);
  for(size_t i = 0; i < n; ++i){
    for(size_t j = i; j < n; ++j){
      int r = rand() % 10000;
      reference_matrix[i][j] = r;
      triangle_matrix.set(i, j, r);
    }
  }
  for(size_t i = 0; i < n; ++i){
    for(size_t j = i+1; j < n; ++j){
      EXPECT_EQ(triangle_matrix.get(i, j), reference_matrix[i][j]);
      EXPECT_EQ(triangle_matrix.get(j, i), reference_matrix[i][j]);
    }
  }
}

TEST_F(TraingleMatrixTest, simple_test){
  size_t n = 10;
  std::vector<std::vector<size_t>> reference_matrix = std::vector<std::vector<size_t>>(n, std::vector<size_t>(n));
  TriangleMatrix<size_t> triangle_matrix = TriangleMatrix<size_t>(n, false);
  for(size_t i = 0; i < n; ++i){
    for(size_t j = i+1; j < n; ++j){
      int r = rand() % 10000;
      reference_matrix[i][j] = r;
      triangle_matrix.set(i, j, r);
    }
  }
  for(size_t i = 0; i < n; ++i){
    for(size_t j = i+1; j < n; ++j){
      EXPECT_EQ(triangle_matrix.get(i, j), reference_matrix[i][j]);
      EXPECT_EQ(triangle_matrix.get(j, i), reference_matrix[i][j]);
    }
  }
}
