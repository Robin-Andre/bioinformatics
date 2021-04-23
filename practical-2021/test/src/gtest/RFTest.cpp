#include "gtest/gtest.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>


class RFTest: public testing::Test{
protected:
  float getDistanceFromString(const std::string &line) const;
  std::vector<float> readDistances(std::string data_set_name) const;
  float getDistance(u_int64_t tree1, u_int64_t tree2, const std::vector<float> &distances, u_int64_t tree_count) const;


};

float RFTest::getDistanceFromString(const std::string &line) const {
    std::istringstream iss (line);
    std::string item;
    u_int64_t i = 0;
    while (std::getline(iss, item, ' ') && i < 3) {
        i++;
    }
    return std::stof(item);
}

 std::vector<float> RFTest::readDistances(std::string data_set_name) const
{
    std::vector<float> distances;
    std::fstream res_file;
    res_file.open("../../../../test/res/reference_results/" + data_set_name + "/RAxML_RF-Distances.0"  ,std::ios::in);
    if (res_file.is_open()){
      std::string line;
      std::vector<std::string> parts;
      while(std::getline(res_file, line)){
        distances.push_back(getDistanceFromString(line));
      }
      res_file.close(); //close the file object.
  }
  return distances;
}

float RFTest::getDistance(u_int64_t tree1, u_int64_t tree2, const std::vector<float> &distances, u_int64_t tree_count) const
{
  if (tree1 == tree2){
    return 0;
  }
  u_int64_t first_tree = std::min(tree1, tree2);
  u_int64_t second_tree = std::max(tree1, tree2);
  u_int64_t offset =(first_tree*(2*tree_count-first_tree-1))/2;
  u_int64_t index = offset + (second_tree - first_tree - 1);
  return distances[index];

}


TEST_F(RFTest, basic_test)
{
    u_int64_t tree_count = 10;
    std::vector<float> distances = readDistances("24");
    std::cout << getDistance(0,2,distances,tree_count) << std::endl;
    EXPECT_EQ(1, 1);

}
