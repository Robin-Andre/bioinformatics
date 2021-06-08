#include "gtest/gtest.h"
#include "../../../src/GeneralizedRFDistance.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/io/IOData.hpp"

class RFTest : public testing::Test {
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
}
void run_test(const std::string& test_file) {
  io::IOData result = GeneralizedRFDistance::computeDistances(
                      TreeReader::readTreeFile(current_data_dir + test_file), rf, ABSOLUTE);
  io::IOData reference = RAXMLReader::read(current_ref_dir + test_file);
  //TODO: Return of Reader is inconstistent - Tip count needs to be set and conversion done in reader
  reference.mean_dst *= (2*(PllSplit::getTipCount()-3));
    ASSERT_EQ(result, reference);
}
};

/*
Yes it would be smart to simply loop over an array of strings and call them instead of this repetitive nonsense
but I would like to have all tests separate
*/
TEST_F(RFTest, 24taxa) {
    execute_test("heads/24");
}
TEST_F(RFTest, 125taxa) {
    execute_test("heads/125");
}
TEST_F(RFTest, 141taxa) {
    execute_test("heads/141");
}
/*TEST_F(RFTest, 143taxa) {
    execute_test("heads/143");
}
TEST_F(RFTest, 148taxa) {
    execute_test("heads/148");
}
TEST_F(RFTest, 150taxa) {
    execute_test("heads/150");
}
TEST_F(RFTest, 218taxa) {
    execute_test("heads/218");
}
TEST_F(RFTest, 350taxa) {
    execute_test("heads/350");
}
TEST_F(RFTest, 354taxa) {
    execute_test("heads/354");
}
TEST_F(RFTest, 404taxa) {
    execute_test("heads/404");
}
TEST_F(RFTest, 500taxa) {
    execute_test("heads/500");
}
TEST_F(RFTest, 628taxa) {
    execute_test("heads/628");
}
TEST_F(RFTest, 714taxa) {
    execute_test("heads/24");
}
TEST_F(RFTest, 885taxa) {
    execute_test("heads/885");
}
TEST_F(RFTest, 994taxa) {
    execute_test("heads/994");
}
TEST_F(RFTest, 1288taxa) {
    execute_test("heads/1288");
}
TEST_F(RFTest, 1481taxa) {
    execute_test("heads/1481");
}
TEST_F(RFTest, 1512taxa) {
    execute_test("heads/1512");
}
TEST_F(RFTest, 1604taxa) {
    execute_test("heads/1604");
}
TEST_F(RFTest, 1908taxa) {
    execute_test("heads/1908");
}
TEST_F(RFTest, 2000taxa) {
    execute_test("heads/2000");
}
TEST_F(RFTest, 2308taxa) {
    execute_test("heads/2308");
}
TEST_F(RFTest, 2554taxa) {
    execute_test("heads/2554");
}*/
