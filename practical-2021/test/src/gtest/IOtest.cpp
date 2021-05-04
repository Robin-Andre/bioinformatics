#include "gtest/gtest.h"
#include "../../../src/io/FileReader.hpp"
#include "../../../src/RFDistance2.hpp"
#include <bitset>
class IOTest : public testing::Test {

};

TEST_F(IOTest, readin) {
  std::vector<PllTree> current = io::readTreeFile("../test/res/data/simple_newick");
  //EXPECT_EQ(current.size(), 10);
  PllSplitList murks = current[0].makeSplits();
  PllSplit test = murks[0];
  PllSplit test2 = murks[1];
  EXPECT_FALSE(test==test2);
  EXPECT_TRUE(test < test2);
}
