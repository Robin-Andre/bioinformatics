#include "gtest/gtest.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <RFDistance.hpp>


class RFTest: public testing::Test{
protected:
  float getDistanceFromString(const std::string &line) const;
  unsigned int readTreeCount(std::string data_set_name) const;
  unsigned int readUniqueTreeCount(std::string data_set_name) const;
  float readAverageDistance(std::string data_set_name) const;
  std::string readFromInfoFile(std::string data_set_name, std::string prefix) const;
  std::vector<float> readDistances(std::string data_set_name) const;
  float getDistance(unsigned int tree1, unsigned int tree2, const std::vector<float> &distances, unsigned int tree_count) const;


};

float RFTest::getDistanceFromString(const std::string &line) const {
    std::istringstream iss (line);
    std::string item;
    unsigned int i = 0;
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

unsigned int RFTest::readTreeCount(std::string data_set_name) const
{
  return std::stoi(readFromInfoFile(data_set_name, "Found "));
}

unsigned int RFTest::readUniqueTreeCount(std::string data_set_name) const
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
   res_file.open("../../../../test/res/reference_results/" + data_set_name + "/RAxML_info.0"  ,std::ios::in);
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

float RFTest::getDistance(unsigned int tree1, unsigned int tree2, const std::vector<float> &distances, unsigned int tree_count) const
{
  if (tree1 == tree2){
    return 0;
  }
  unsigned int first_tree = std::min(tree1, tree2);
  unsigned int second_tree = std::max(tree1, tree2);
  unsigned int offset =(first_tree*(2*tree_count-first_tree-1))/2;
  unsigned int index = offset + (second_tree - first_tree - 1);
  return distances[index];

}


TEST_F(RFTest, basic_test)
{
    std::string test_set = "24";
    float error = 0.01;

    RFDistance rf_distance = RFDistance("../../../../test/res/data/heads/BS/" + test_set);
    rf_distance.run();
    rf_distance.writeResults("../../../../output/" + test_set);
    EXPECT_EQ(rf_distance.getTreeCount(), readTreeCount(test_set));
    unsigned int tree_count = rf_distance.getTreeCount();
    EXPECT_EQ(rf_distance.getUniqueCount(), readTreeCount(test_set));
    EXPECT_NEAR(rf_distance.getAverageDistance(), readAverageDistance(test_set), error);
    unsigned int k=0;
    std::vector<float> calculated_distances = rf_distance.getDistances();
    std::vector<float> reference_distances = readDistances(test_set);
    for (unsigned int i=0; i < tree_count; i++){
      for (unsigned int j=i+1; j < tree_count; j++){
        std::cout << k << std::endl;
        EXPECT_NEAR(calculated_distances[k], reference_distances[k], error);
        k++;
      }
    }

    //std::cout << getDistance(0,2,distances,tree_count) << std::endl;


}
