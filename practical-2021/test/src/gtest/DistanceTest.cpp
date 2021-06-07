#include "gtest/gtest.h"
#include "../../../src/GeneralizedRFDistance.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/io/IOData.hpp"

class DistanceTest : public testing::Test {
protected:
/*This is a hardcoded link to the test dir. IF changes to the project structure are made this needs
to be adjusted.
*/
//std::string current_test_dir = "../test/res/data/heads/BS/";
std::string current_data_dir = "../test/res/data/";
std::string current_ref_dir = "../test/res/reference_results/";
float epsilon = 0.001;
RFMetric rf;

/*Method to reduce code complexity :)
*/
void execute_test(const std::string& test_file) {
    run_test(test_file);
    ensure_no_bitmask_needed(test_file);
}
void run_test(const std::string& test_file) {
  io::IOData result = GeneralizedRFDistance::computeDistances(
                      TreeReader::readTreeFile(current_data_dir + test_file), rf, false);
  io::IOData reference = RAXMLReader::read(current_ref_dir + test_file);
    ASSERT_EQ(result, reference);
}
/*I am still convinced that the bitmask of the PllSplits is weird so I checked that no PllSplit
   requires that mask.
  The original implementation was bitmask on register 0. Should have been register n.
  Ok turns out that the PllSplitList object needs to be held to have a reasonable splitlist.
*/
void ensure_no_bitmask_needed(const std::string& test_file) {

    PllTree test_tree = TreeReader::readTreeFile(current_data_dir + test_file)[0];
    //this shit works
    PllSplitList object = PllSplitList(test_tree);
    std::vector<PllSplit> splits = object.getSplits();

    //this shit is broken
    //std::vector<PllSplit> splits = PllSplitList(test_tree).getSplits();

    PllSplit::setTipCount(test_tree.getTipCount());
    for(size_t i = 0; i < splits.size(); ++i) {
      ASSERT_EQ(splits[i].bitExtract(0), 1); //First bit is always 1
      for(size_t j = test_tree.getTipCount(); j < PllSplit::getSplitLen() * sizeof(pll_split_base_t) * 8; ++j) {
        ASSERT_EQ(splits[i].bitExtract(j), 0);
      }
    }


}
};
/*
Yes it would be smart to simply loop over an array of strings and call them instead of this repetitive nonsense
but I would like to have all tests separate
*/
TEST_F(DistanceTest, 24taxa) {
    execute_test("heads/24");
}
TEST_F(DistanceTest, 125taxa) {
    execute_test("heads/125");
}
TEST_F(DistanceTest, 141taxa) {
    execute_test("heads/141");
}
/*TEST_F(DistanceTest, 143taxa) {
    execute_test("heads/143");
}
TEST_F(DistanceTest, 148taxa) {
    execute_test("heads/148");
}
TEST_F(DistanceTest, 150taxa) {
    execute_test("heads/150");
}
TEST_F(DistanceTest, 218taxa) {
    execute_test("heads/218");
}
TEST_F(DistanceTest, 350taxa) {
    execute_test("heads/350");
}
TEST_F(DistanceTest, 354taxa) {
    execute_test("heads/354");
}
TEST_F(DistanceTest, 404taxa) {
    execute_test("heads/404");
}
TEST_F(DistanceTest, 500taxa) {
    execute_test("heads/500");
}
TEST_F(DistanceTest, 628taxa) {
    execute_test("heads/628");
}
TEST_F(DistanceTest, 714taxa) {
    execute_test("heads/24");
}
TEST_F(DistanceTest, 885taxa) {
    execute_test("heads/885");
}
TEST_F(DistanceTest, 994taxa) {
    execute_test("heads/994");
}
TEST_F(DistanceTest, 1288taxa) {
    execute_test("heads/1288");
}
TEST_F(DistanceTest, 1481taxa) {
    execute_test("heads/1481");
}
TEST_F(DistanceTest, 1512taxa) {
    execute_test("heads/1512");
}
TEST_F(DistanceTest, 1604taxa) {
    execute_test("heads/1604");
}
TEST_F(DistanceTest, 1908taxa) {
    execute_test("heads/1908");
}
TEST_F(DistanceTest, 2000taxa) {
    execute_test("heads/2000");
}
TEST_F(DistanceTest, 2308taxa) {
    execute_test("heads/2308");
}
TEST_F(DistanceTest, 2554taxa) {
    execute_test("heads/2554");
}*/
