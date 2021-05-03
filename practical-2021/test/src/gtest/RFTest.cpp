#include "gtest/gtest.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <../../../src/RFDistance.hpp>


class RFTest: public testing::Test{
protected:
  size_t getDistanceFromString(const std::string &line) const;
  size_t readTreeCount(std::string data_set_name) const;
  size_t readUniqueTreeCount(std::string data_set_name) const;
  float readAverageDistance(std::string data_set_name) const;
  std::string readFromInfoFile(std::string data_set_name, std::string prefix) const;
  std::vector<size_t> readDistances(std::string data_set_name) const;
  size_t getDistance(size_t tree1, size_t tree2, const std::vector<size_t> &distances, size_t tree_count) const;


};

size_t RFTest::getDistanceFromString(const std::string &line) const {
    std::istringstream iss (line);
    std::string item;
    size_t i = 0;
    while (std::getline(iss, item, ' ') && i < 2) {
        i++;

    }
    return std::stoi(item);
}

 std::vector<size_t> RFTest::readDistances(std::string data_set_name) const
{
    std::vector<size_t> distances;
    std::fstream res_file;
    res_file.open("../test/res/reference_results/" + data_set_name + "/RAxML_RF-Distances.0"  ,std::ios::in);
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

size_t RFTest::readTreeCount(std::string data_set_name) const
{
  return std::stoi(readFromInfoFile(data_set_name, "Found "));
}

size_t RFTest::readUniqueTreeCount(std::string data_set_name) const
{
  return std::stoi(readFromInfoFile(data_set_name, "Number of unique trees in this tree set: "));
}


float RFTest::readAverageDistance(std::string data_set_name) const
{
  return std::stof(readFromInfoFile(data_set_name, "Average relative RF in this set: "));
}

std::string RFTest::readFromInfoFile(std::string data_set_name, std::string prefix) const
{
   std::string result;
   std::fstream res_file;
   res_file.open("../test/res/reference_results/" + data_set_name + "/RAxML_info.0"  ,std::ios::in);
   if (res_file.is_open()){
     std::string line;
     std::vector<std::string> parts;
     while(std::getline(res_file, line)){
       auto match = std::mismatch(prefix.begin(), prefix.end(), line.begin());
       if (match.first == prefix.end())
       {
         std::istringstream iss (line);
         std::istringstream iss2 (prefix);
         std::string dummy;
         while (std::getline(iss2, dummy, ' ')){
           std::getline(iss, result, ' ');
         }
         std::getline(iss, result, ' ');
         break;
       }
     }
     res_file.close(); //close the file object.
 }
 return result;
}

size_t RFTest::getDistance(size_t tree1, size_t tree2, const std::vector<size_t> &distances, size_t tree_count) const
{
  if (tree1 == tree2){
    return 0;
  }
  size_t first_tree = std::min(tree1, tree2);
  size_t second_tree = std::max(tree1, tree2);
  size_t offset =(first_tree*(2*tree_count-first_tree-1))/2;
  size_t index = offset + (second_tree - first_tree - 1);
  return distances[index];

}


TEST_F(RFTest, basic_test)
{
    std::string test_set = "350";
    float error = 0.01;
    RFDistance rf_distance = RFDistance();
    rf_distance.run("../test/res/data/heads/BS/" + test_set);
    rf_distance.writeResults("../output/" + test_set);
    EXPECT_EQ(rf_distance.getTreeCount(), readTreeCount(test_set));
    size_t tree_count = rf_distance.getTreeCount();
    EXPECT_EQ(rf_distance.getUniqueCount(), readUniqueTreeCount(test_set));
    EXPECT_NEAR(rf_distance.getAverageDistance(), readAverageDistance(test_set), error);
    size_t k=0;
    std::vector<size_t> calculated_distances = rf_distance.getDistances();
    std::vector<size_t> reference_distances = readDistances(test_set);
    for (size_t i=0; i < tree_count; i++){
      for (size_t j=i+1; j < tree_count; j++){
        EXPECT_EQ(calculated_distances[k], reference_distances[k]) << i << " " << j;
        k++;
      }
    }

    //std::cout << getDistance(0,2,distances,tree_count) << std::endl;


}
