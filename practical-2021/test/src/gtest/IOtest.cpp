#include "gtest/gtest.h"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataWriter.hpp"
#include "../TestUtil.hpp"
#include <filesystem>
class IOTest : public testing::Test {
  float epsilon = 0.00000001;

};

TEST_F(IOTest, tree_read) {
  PllSplitList splits = TreeReader::readTreeFile("../test/res/data/simple_newick")[0];
  PllSplit split1 = splits[0];
  PllSplit split2 = splits[1];
  EXPECT_FALSE(split1==split2);
  EXPECT_TRUE(split1 < split2);
}

TEST_F(IOTest, raxml_read_write) {
  std::filesystem::create_directories("./foo");
  RFData input = RAXMLReader::read("../test/res/reference_results/heads/24");
  RAXMLWriter::write("foo" , input);
  TestUtil::rf_data_eq(input, RAXMLReader::read("foo"), 0);
  std::filesystem::remove_all("./foo");

}
