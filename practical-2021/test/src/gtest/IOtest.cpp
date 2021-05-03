#include "gtest/gtest.h"
#include "../../../src/io/FileReader.hpp"
class IOTest : public testing::Test {

};

TEST_F(IOTest, readin) {
  std::vector<PllTree> current = io::readTreeFile("../test/res/data/heads/BS/24");
  EXPECT_EQ(current.size(), 10);
  PllSplitList murks = current[0].makeSplits();
  std::cout << murks[0].popcount() << "\n";
  EXPECT_EQ(1,1);
}
