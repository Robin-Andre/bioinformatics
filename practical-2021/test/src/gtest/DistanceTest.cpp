#include "gtest/gtest.h"
#include "../../../src/RFDistance2.hpp"
#include "../../../src/RFDistance.hpp"
#include "../../../src/io/FileReader.hpp"
class DistanceTest : public testing::Test {
protected:
/*Right now an instanciation of test is needed, if we turn it into a free function this needs
to be adapted*/
RFDistance test;
/*This is a hardcoded link to the test dir. IF changes to the project structure are made this needs
to be adjusted.
*/
//std::string current_test_dir = "../test/res/data/heads/BS/";
std::string current_test_dir = "../test/res/data/";
float epsilon = 0.001;
/*
  Compares two vectors, there is probably a smart method to do this within googletest
  but I don't know it (yet). Maybe a macro, maybe a selfdefined template...
*/
void evaluate(const RFData& result, const std::vector<size_t>& expected_values) {
    EXPECT_EQ(result.distances.size(), expected_values.size());
    for(unsigned i = 0; i < expected_values.size(); ++i) {
        EXPECT_EQ(result.distances[i], expected_values[i]);
    }
}
/*Method to reduce code complexity :)
*/
void execute_test(std::string test_file) {
    RFData results = test.computeRF(current_test_dir + test_file);
    evaluate(results, io::readDistances(test_file));
    EXPECT_NEAR(results.average_distance, io::readAverageDistance(test_file), epsilon);
    EXPECT_EQ(results.unique_count, io::readUniqueTreeCount(test_file));

}
};
/*
Yes it would be smart to simply loop over an array of strings and call them instead of this repetitive nonsense
but I would like to have all tests separate
*/
TEST_F(DistanceTest, 24taxa) {
    execute_test("heads/24");
<<<<<<< HEAD
=======
}
/*TEST_F(DistanceTest, 24taxa) {
    execute_test("heads/24");
>>>>>>> d41fbb94b4c4a1b744b17407c0a08370f5111400
}
/*TEST_F(DistanceTest, 125taxa) {
    execute_test("heads/125");
}
TEST_F(DistanceTest, 141taxa) {
    execute_test("heads/141");
}
TEST_F(DistanceTest, 143taxa) {
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
