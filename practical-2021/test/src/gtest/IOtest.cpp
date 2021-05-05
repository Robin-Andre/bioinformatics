#include "gtest/gtest.h"
#include "../../../src/io/FileReader.hpp"
#include <bitset>
class IOTest : public testing::Test {

};

TEST_F(IOTest, readin) {
  PllSplitList splits = io::readTreeFile("../test/res/data/simple_newick")[0];
  PllSplit split1 = splits[0];
  PllSplit split2 = splits[1];
  EXPECT_FALSE(split1==split2);
  EXPECT_TRUE(split1 < split2);
}
