#include "gtest/gtest.h"
#include "io/TreeReader.hpp"
#include "io/RFDataReader.hpp"
#include "io/RFDataWriter.hpp"
#include "io/IOData.hpp"
#include <filesystem>
class IOTest : public testing::Test {
  float epsilon = 0.00000001;

};

/**
* Test tree reader with a sample file
*/

TEST_F(IOTest, tree_read) {
  PllTree tree = TreeReader::readTreeFile("../test/res/data/simple_newick")[0];
  PllSplit::setTipCount(tree.getTipCount());
  UniquePllMap map({tree});
  PllSplit split1 = map[0];
  PllSplit split2 = map[1];
  EXPECT_FALSE(split1==split2);
  EXPECT_TRUE(split1 < split2);
}

/**
*
* Test reading and writing JSOM
*/
TEST_F(IOTest, json_read_write) {
  std::filesystem::create_directories("./foo");
  io::IOData input = JSONReader::read("../test/res/references_json/RF/ABSOLUTE/heads/24");
  JSONWriter::write("foo/24" , input);
  EXPECT_EQ(input, JSONReader::read("foo/24"));
  std::filesystem::remove_all("./foo");
}

//Files of this format are no longer in repo, code ist kept for conversion
/*TEST_F(IOTest, raxml_read_write) {
  std::filesystem::create_directories("./foo");
  io::IOData input = RAXMLReader::read("../test/res/reference_results/heads/24");
  RAXMLWriter::write("foo" , input);
  EXPECT_EQ(input, RAXMLReader::read("foo"));
  std::filesystem::remove_all("./foo");
}

TEST_F(IOTest, matrix_read_write) {
  std::filesystem::create_directories("./foo");
  io::IOData input = MatrixReader::read("../test/res/R_results/MCI/SIMILARITY/heads/24");
  MatrixWriter::write("foo/24" , input);
  EXPECT_EQ(input, MatrixReader::read("foo/24"));
  std::filesystem::remove_all("./foo");
}*/
