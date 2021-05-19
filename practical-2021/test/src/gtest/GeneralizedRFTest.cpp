#include "gtest/gtest.h"
#include "../../../src/DistanceUtil.hpp"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataReader.hpp"

class GeneralizedRFTest : public testing::Test {
protected:
/*Right now an instanciation of test is needed, if we turn it into a free function this needs
to be adapted*/
//GeneralizedRFTest test;
/*This is a hardcoded link to the test dir. IF changes to the project structure are made this needs
to be adjusted.
*/
//std::string current_test_dir = "../test/res/data/heads/BS/";
std::string current_data_dir = "../test/res/data/";
std::string current_ref_dir = "../test/res/reference_results/";
float epsilon = 0.001;

/*Method to reduce code complexity :)
*/
void execute_test(std::string test_file) {
    std::vector<PllTree> trees = TreeReader::readTreeFile(current_data_dir + test_file);
    PllSplit::setTipCount(trees[0].getTipCount());
    PllSplitList first = PllSplitList(trees[0]);
    PllSplitList second = PllSplitList(trees[1]);
    std::vector<std::vector<double>> distances = DistanceUtil::similaritiesForSplits(first, second, SPI);
    for(size_t i = 0; i < distances.size(); ++i){
      for(size_t j = 0; j < distances[i].size(); ++j){
        std::cout << distances[i][j] << ";";
      }
      std::cout << std::endl;
    }
}
};
/*
Yes it would be smart to simply loop over an array of strings and call them instead of this repetitive nonsense
but I would like to have all tests separate
*/
TEST_F(GeneralizedRFTest, 24taxa) {
    execute_test("heads/24");
}
/*TEST_F(GeneralizedRFTest, 125taxa) {
    execute_test("heads/125");
}
TEST_F(GeneralizedRFTest, 141taxa) {
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
