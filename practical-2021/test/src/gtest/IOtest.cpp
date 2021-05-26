#include "gtest/gtest.h"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataWriter.hpp"
#include <filesystem>
class IOTest : public testing::Test {
  float epsilon = 0.00000001;

};

TEST_F(IOTest, tree_read) {
  PllTree tree = TreeReader::readTreeFile("../test/res/data/simple_newick")[0];
  PllSplit::setTipCount(tree.getTipCount());
  PllSplitList splits = PllSplitList(tree);
  PllSplit split1 = splits[0];
  PllSplit split2 = splits[1];
  EXPECT_FALSE(split1==split2);
  EXPECT_TRUE(split1 < split2);
}

TEST_F(IOTest, raxml_read_write) {
  std::filesystem::create_directories("./foo");
  RFData input = RAXMLReader::read("../test/res/reference_results/heads/24");
  RAXMLWriter::write("foo" , input);
  EXPECT_EQ(input, RAXMLReader::read("foo"));
  std::filesystem::remove_all("./foo");
}

TEST_F(IOTest, json_read_write) {
  std::filesystem::create_directories("./foo");
  RFData input = RAXMLReader::read("../test/res/reference_results/heads/24");
  JSONWriter::write("foo" , input);
  EXPECT_EQ(input, JSONReader::read("foo"));
  std::filesystem::remove_all("./foo");
}
/*TEST_F(IOTest, proper_print) {
  PllTree tree = TreeReader::readTreeFile("../test/res/data/heads/125")[0];
  PllSplit::setTipCount(tree.getTipCount());
  PllSplitList splits = PllSplitList(tree);
  splits.printSplits();
}*/
TEST_F(IOTest, is_the_number_of_tips_accurate) {
  PllTree tree = TreeReader::readTreeFile("../test/res/data/heads/125")[0];
  ASSERT_EQ(tree.getTipCount(), 125);
}