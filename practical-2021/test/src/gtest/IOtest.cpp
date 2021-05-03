#include "gtest/gtest.h"
#include "../../../src/io/FileReader.hpp"
#include "../../../src/RFDistance2.hpp"
#include <bitset>
class IOTest : public testing::Test {

};

TEST_F(IOTest, readin) {
  std::vector<PllTree> current = io::readTreeFile("../test/res/data/heads/BS/simple_newick");
  //EXPECT_EQ(current.size(), 10);
  PllSplitList murks = current[0].makeSplits();
  PllSplit test = murks[0];
  PllSplit test2 = murks[1];
  std::cout << murks[0].popcount() << "\n";
  std::cout << (test == test2) << (test < test2) << "EQUALITY\n";
  std::cout << *test() << "\n";
  for(unsigned i = 0; i < 24; ++i) {
      std::cout << "Popcount with: " << i << " elements: " << test.popcount() << "\n";
  }
  pll_split_base_t lol = test()[0];
  std::cout << lol << "OFLF\n";
  std::cout << sizeof(pll_split_base_t) * 8 << " THATSTHESIZE\n";
   std::cout << *(test()+1) << "\n";
   std::cout << *(test()+2) << "\n";
  EXPECT_EQ(1,1);
}
