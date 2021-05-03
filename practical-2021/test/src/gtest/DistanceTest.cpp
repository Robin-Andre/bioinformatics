#include "gtest/gtest.h"
#include "../../../src/RFDistance2.hpp"
#include "../../../src/RFDistance.hpp"
#include "../../../src/io/FileReader.hpp"
class DistanceTest : public testing::Test {
protected:
RFDistance test;
std::string current_test_dir = "../test/res/data/heads/BS/";
void evaluate(const std::vector<size_t>& actual_values, const std::vector<size_t>& expected_values) {
    EXPECT_EQ(actual_values.size(), expected_values.size());
    for(unsigned i = 0; i < actual_values.size(); ++i) {
        EXPECT_EQ(actual_values[i], expected_values[i]);
    }
}
};

TEST_F(DistanceTest, 24taxa) {
    std::string test_file = "24";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 125taxa) {
    std::string test_file = "125";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 141taxa) {
    std::string test_file = "141";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 143taxa) {
    std::string test_file = "143";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 148taxa) {
    std::string test_file = "148";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 150taxa) {
    std::string test_file = "150";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 218taxa) {
    std::string test_file = "218";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 350taxa) {
    std::string test_file = "350";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 354taxa) {
    std::string test_file = "354";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 404taxa) {
    std::string test_file = "404";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 500taxa) {
    std::string test_file = "500";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 628taxa) {
    std::string test_file = "628";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 714taxa) {
    std::string test_file = "714";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 885taxa) {
    std::string test_file = "885";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 994taxa) {
    std::string test_file = "994";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 1288taxa) {
    std::string test_file = "1288";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 1481taxa) {
    std::string test_file = "1481";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 1512taxa) {
    std::string test_file = "1512";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 1604taxa) {
    std::string test_file = "1604";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 1908taxa) {
    std::string test_file = "1908";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 2000taxa) {
    std::string test_file = "2000";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 2308taxa) {
    std::string test_file = "2308";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
TEST_F(DistanceTest, 2554taxa) {
    std::string test_file = "2554";
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results.distances, io::readDistances(test_file));
}
