#include "gtest/gtest.h"
#include "../../../src/io/FileReader.hpp"
#include "../../../src/RFDistance2.hpp"
#include <bitset>
class IOTest : public testing::Test {

};

<<<<<<< HEAD
TEST_F(IOTest, readin) {
  std::vector<PllTree> current = io::readTreeFile("../test/res/data/simple_newick");
  //EXPECT_EQ(current.size(), 10);
  PllSplitList murks = current[0].makeSplits();
  PllSplit test = murks[0];
  PllSplit test2 = murks[1];
  EXPECT_FALSE(test==test2);
  EXPECT_TRUE(test < test2);
}
=======
/*TEST_F(IOTest, readin) {
  PllTree simple_tree = io::readTreeFile("../test/res/data/simple_newick")[0];
  PllSplitList simple_splits = simple_tree.makeSplits();
  PllSplit test = simple_splits[0];
  PllSplit test2 = simple_splits[1];
  std::cout << "Popcount of first split "<<simple_splits[0].popcount() << "\n";
  std::cout << (simple_splits[0] == simple_splits[1]) << (simple_splits[0] < simple_splits[1]) << " EQUALITY\n";
  /*std::cout << *test() << "\n";
  for(unsigned i = 0; i < 24; ++i) {
      std::cout << "Popcount with: " << i << " elements: " << test.popcount() << "\n";
  }
   std::cout << *(test()+1) << "\n";
   std::cout << *(test()+2) << "\n";
}*/
>>>>>>> d41fbb94b4c4a1b744b17407c0a08370f5111400
