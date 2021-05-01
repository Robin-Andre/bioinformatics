#include "gtest/gtest.h"
#include "../../../src/io/FileReader.hpp"
#include "../../../src/RFDistance2.hpp"
#include <bitset>
class IOTest : public testing::Test {
    
};
TEST_F(IOTest, full_calculation) {
    //The results are actually correct 
    std::vector<size_t> results = measure::full_calculation("../test/res/data/heads/BS/24");
    std::cout << results.size() << "VALUES THAT I DON'T CARE ABOUT\n";

}
TEST_F(IOTest, readin) {
  std::vector<PllTree> current = io::readTreeFile("../test/res/data/heads/BS/simple_newick");
  //EXPECT_EQ(current.size(), 10);
  PllSplitList murks = current[0].makeSplits();
  
  
  PllSplit test = murks[0];
  PllSplit test2 = murks[1];
  std::cout << test.popcount_vector() << "EGEGEGEG\n";
  std::cout << murks[0].popcount(1) << "\n";
  std::cout << (test == test2) << (test < test2) << "EQUALITY\n";
  std::cout << std::bitset<5>(*test()) <<  " " <<std::bitset<5>(*test2()) << "\n";
  std::cout << test.bitExtract(0) << " "<< test.bitExtract(1) << test.bitExtract(2) << test.bitExtract(3) << test.bitExtract(4)  << " " << "WOLOLOL" << "\n";
  std::cout << test2.bitExtract(0) << " "<< test2.bitExtract(1) << test2.bitExtract(2) << test2.bitExtract(3) << test2.bitExtract(4)  << " " << "WOLOLOL" << "\n";
  std::cout << *test() << "\n";
  for(unsigned i = 0; i < 24; ++i) {
      std::cout << "Popcount with: " << i << " elements: " << test.popcount(i) << "\n";
  }
  pll_split_base_t lol = test()[0];
  std::cout << lol << "OFLF\n";
  std::cout << sizeof(pll_split_base_t) * 8 << "THATSTHESIZE\n";
  std::cout << test.popcount(1) << "\n";
   std::cout << *(test()+1) << "\n";
   std::cout << *(test()+2) << "\n";
  EXPECT_EQ(1,1);
}
// These are from the 24 set, tree 0 and 2 respectively (Alexis code says distance 4)
TEST_F(IOTest, UnequalTrees) {
    std::vector<PllTree> current = io::readTreeFile("../test/res/data/heads/BS/TwoTreesUnequal");
    current[0].alignNodeIndices(current[1]);
    PllSplitList firstset = current[0].makeSplits();
    PllSplitList secondset = current[1].makeSplits();
    
    std::cout << measure::rf_distance(firstset, secondset) << "Thisisaresult and it might be shit\n";
    EXPECT_EQ(measure::rf_distance(firstset, secondset) ,4);
}
TEST_F(IOTest, EqualTrees) {
    std::vector<PllTree> current = io::readTreeFile("../test/res/data/heads/BS/TwoTreesEqual");
    PllSplitList firstset = current[0].makeSplits();
    PllSplitList secondset = current[1].makeSplits();
    current[0].alignNodeIndices(current[1]);
    //std::cout << measure::rf_distance(firstset, secondset) << "Thisisaresult and it might be shit\n";
    EXPECT_EQ(measure::rf_distance(firstset, secondset),0);
}

//Remember, the lecture slide is wrong and the RF distance is actually 2. 
TEST_F(IOTest, LectureExample) {
    std::vector<PllTree> current = io::readTreeFile("../test/res/data/heads/BS/example_from_slideshow");
    current[0].alignNodeIndices(current[1]);
    PllSplitList firstset = current[0].makeSplits();
    PllSplitList secondset = current[1].makeSplits();
    std::cout << *firstset[0]() << " " << *firstset[1]() << " " << *firstset[2]() << "Thefirstset\n";
    std::cout << std::bitset<6>(*firstset[0]())<< " " << std::bitset<6>(*firstset[1]()) << " " << std::bitset<6>(*firstset[2]()) << "Thefirstset\n";
    std::cout << std::bitset<6>(*secondset[0]())<< " " << std::bitset<6>(*secondset[1]()) << " " << std::bitset<6>(*secondset[2]()) << "Thefirstset\n";
    //std::cout << measure::rf_distance(firstset, secondset) << "Thisisaresult and it might be shit\n";
    EXPECT_EQ(measure::rf_distance(firstset, secondset),2);
}