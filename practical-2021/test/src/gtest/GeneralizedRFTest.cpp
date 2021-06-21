#include "gtest/gtest.h"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataWriter.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/io/IOData.hpp"
#include "../../../src/GeneralizedRFDistance.hpp"
#include "../../../src/Metric.hpp"

#include <ortools/linear_solver/linear_solver.h>
#include <string>

static constexpr bool print_execution_time = false;

using GRFDist = GeneralizedRFDistance;
class GeneralizedRFTest : public testing::Test {

protected:
/*Right now an instanciation of test is needed, if we turn it into a free function this needs
to be adapted*/
std::vector<PllTree> load_24taxa() {
   PllSplit::setTipCount(24);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/24")[1];
  std::vector<PllTree> trees = {tree1, tree2};
  return trees;
}
//GeneralizedRFTest test;
/*This is a hardcoded link to the test dir. IF changes to the project structure are made this needs
to be adjusted.
*/
//std::string current_test_dir = "../test/res/data/heads/BS/";
std::string current_data_dir = "../test/res/data/";
std::string current_ref_dir = "../test/res/references_json/";
float epsilon = 0.001;

SPIMetric metric_spi;
MCIMetric metric_mci;
MSIMetric metric_msi;
RFMetric metric_rf;

/*Method to reduce code complexity :)
*/
void execute_test(const std::string& test_file, const Metric& metric, Mode mode) {
    std::string mode_name = ModeString[mode];
    std::vector<PllTree> trees = TreeReader::readTreeFile(current_data_dir + test_file);
    auto time_start = std::chrono::high_resolution_clock::now();
    io::IOData result = GeneralizedRFDistance::computeDistances(trees, metric, mode);
    auto time_end = std::chrono::high_resolution_clock::now();
    io::IOData reference = JSONReader::read(current_ref_dir + metric.name() + "/" + mode_name + "/" + test_file);
    EXPECT_EQ(result, reference);
    //This is a dirty solution to get execution time on the console for measurement purposes
    //Right now all tests are intertwined and this is the easiest insertion point 
    if(print_execution_time) {
      std::cout << "TESTOUTPUT: " << test_file << " " << metric.name() << " " << mode 
                << " " << (std::chrono::duration<double, std::milli>(time_end - time_start)).count() << "\n";
    }

}


};

TEST_F(GeneralizedRFTest, simple_identity) {
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  std::vector<PllTree> trees = {tree, tree};
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_msi, RELATIVE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_spi, RELATIVE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_mci, RELATIVE).pairwise_distance_mtx[0][0], 0, epsilon);
}
TEST_F(GeneralizedRFTest, ExampleFromSlideshow) {
  PllSplit::setTipCount(6);
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  std::vector<PllTree> trees = {tree, tree};
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}
TEST_F(GeneralizedRFTest, ComparisionTree0_2taxa24) {
  std::vector<PllTree> trees = load_24taxa();
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}
TEST_F(GeneralizedRFTest, example_24_msi) {
  std::vector<PllTree> trees = load_24taxa();
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}

TEST_F(GeneralizedRFTest, example_from_slideshow) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  std::vector<PllTree> trees = {tree1, tree2};
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeDistances(trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}
TEST_F(GeneralizedRFTest, 24taxa) {
  execute_test("heads/24", metric_rf, ABSOLUTE);
  execute_test("heads/24", metric_rf, RELATIVE);

  execute_test("heads/24", metric_msi, SIMILARITY);
  execute_test("heads/24", metric_spi, SIMILARITY);
  execute_test("heads/24", metric_mci, SIMILARITY);
  execute_test("heads/24", metric_mci, ABSOLUTE);
}
/*TEST_F(GeneralizedRFTest, 125taxa) {
  execute_test("heads/125", metric_rf, ABSOLUTE);
  execute_test("heads/125", metric_rf, RELATIVE);

  execute_test("heads/125", metric_msi, SIMILARITY);
  execute_test("heads/125", metric_spi, SIMILARITY);
  execute_test("heads/125", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 141taxa) {
  execute_test("heads/141", metric_msi, SIMILARITY);
  execute_test("heads/141", metric_spi, SIMILARITY);
  execute_test("heads/141", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 143taxa) {
  execute_test("heads/143", metric_msi, SIMILARITY);
  execute_test("heads/143", metric_spi, SIMILARITY);
  execute_test("heads/143", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 148taxa) {
  execute_test("heads/148", metric_msi, SIMILARITY);
  execute_test("heads/148", metric_spi, SIMILARITY);
  execute_test("heads/148", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 150taxa) {
  execute_test("heads/150", metric_msi, SIMILARITY);
  execute_test("heads/150", metric_spi, SIMILARITY);
  execute_test("heads/150", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 218taxa) {
  execute_test("heads/218", metric_msi, SIMILARITY);
  execute_test("heads/218", metric_spi, SIMILARITY);
  execute_test("heads/218", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 350taxa) {
  execute_test("heads/350", metric_msi, SIMILARITY);
  execute_test("heads/350", metric_spi, SIMILARITY);
  execute_test("heads/350", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 354taxa) {
  execute_test("heads/354", metric_msi, SIMILARITY);
  execute_test("heads/354", metric_spi, SIMILARITY);
  execute_test("heads/354", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 404taxa) {
  execute_test("heads/404", metric_msi, SIMILARITY);
  execute_test("heads/404", metric_spi, SIMILARITY);
  execute_test("heads/404", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 500taxa) {
  execute_test("heads/500", metric_msi, SIMILARITY);
  execute_test("heads/500", metric_spi, SIMILARITY);
  execute_test("heads/500", metric_mci, SIMILARITY);
}

TEST_F(GeneralizedRFTest, 885taxa) {
  execute_test("heads/885", metric_msi, SIMILARITY);
  execute_test("heads/885", metric_spi, SIMILARITY);
  execute_test("heads/885", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 994taxa) {
  execute_test("heads/994", metric_msi, SIMILARITY);
  execute_test("heads/994", metric_spi, SIMILARITY);
  execute_test("heads/994", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 1288taxa) {
  execute_test("heads/1288", metric_msi, SIMILARITY);
  execute_test("heads/1288", metric_spi, SIMILARITY);
  execute_test("heads/1288", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 1481taxa) {
  execute_test("heads/1481", metric_msi, SIMILARITY);
  execute_test("heads/1481", metric_spi, SIMILARITY);
  execute_test("heads/1481", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 1512taxa) {
  execute_test("heads/1512", metric_msi, SIMILARITY);
  execute_test("heads/1512", metric_spi, SIMILARITY);
  execute_test("heads/1512", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 1604taxa) {
  execute_test("heads/1604", metric_msi, SIMILARITY);
  execute_test("heads/1604", metric_spi, SIMILARITY);
  execute_test("heads/1604", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 1908taxa) {
  execute_test("heads/1908", metric_msi, SIMILARITY);
  execute_test("heads/1908", metric_spi, SIMILARITY);
  execute_test("heads/1908", metric_mci, SIMILARITY);
}
TEST_F(GeneralizedRFTest, 2000taxa) {
  execute_test("heads/2000", metric_msi, SIMILARITY);
  execute_test("heads/2000", metric_spi, SIMILARITY);
  execute_test("heads/2000", metric_mci, SIMILARITY);
}*/
