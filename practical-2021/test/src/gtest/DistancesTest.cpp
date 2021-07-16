#include "gtest/gtest.h"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataWriter.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/io/IOData.hpp"
#include "../../../src/Distances.hpp"
#include "../../../src/Metric.hpp"

#include <ortools/linear_solver/linear_solver.h>
#include <string>

#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif

static constexpr bool print_execution_time = true;

using GRFDist = Distances;
class DistancesTest : public testing::Test {

protected:
//TODO this should be compile macroed
/*static void SetUpTestSuite() {
  if(print_execution_time) {
    io::clear_benchmark_timing();
  }

}*/
/*Right now an instanciation of test is needed, if we turn it into a free function this needs
to be adapted*/
std::vector<PllTree> load_24taxa() {
   PllSplit::setTipCount(24);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "heads/24")[1];
  std::vector<PllTree> trees = {tree1, tree2};
  return trees;
}
//DistancesTest test;
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

void execute_test_rf(const std::string& test_file, Mode mode) {
  io::IOData result = Distances::computeRFDistances(
                      TreeReader::readTreeFile(current_data_dir + test_file), mode);
  io::IOData reference = JSONReader::read(current_ref_dir + "RF/" + ModeString[mode] + "/" +test_file);
    ASSERT_EQ(result, reference);
}


void execute_test_generalized(const std::string& test_file, const GeneralizedMetric& metric, Mode mode) {
    std::string mode_name = ModeString[mode];
    std::vector<PllTree> trees = TreeReader::readTreeFile(current_data_dir + test_file);
    auto time_start = std::chrono::high_resolution_clock::now();
    io::IOData result = Distances::computeGeneralizedDistances(trees, metric, mode);

    auto time_end = std::chrono::high_resolution_clock::now();
    io::IOData reference = JSONReader::read(
        current_ref_dir + metric.name() + "/" + mode_name + "/" + test_file);
    EXPECT_EQ(result, reference);
    //This is a dirty solution to get execution time on the console for measurement purposes
    //Right now all tests are intertwined and this is the easiest insertion point
    if(print_execution_time) {

      std::string result = test_file + " " + metric.name() + " " + mode_name
                + " " + std::to_string((std::chrono::duration<double, std::milli>(time_end - time_start)).count())
                + " "+ GIT_COMMIT_HASH;
      io::write_benchmark_timing(result);
    }

}


};

TEST_F(DistancesTest, simple_identity) {
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  std::vector<PllTree> trees = {tree, tree};
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_msi, RELATIVE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_spi, RELATIVE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_mci, RELATIVE).pairwise_distance_mtx[0][0], 0, epsilon);

}

TEST_F(DistancesTest, no_unique_tree) {
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "heads/24")[0];
  std::vector<PllTree> trees = {tree, tree, tree, tree};

  EXPECT_EQ(GRFDist::computeGeneralizedDistances(trees, metric_msi, ABSOLUTE).number_of_unique_trees, 0);
  EXPECT_EQ(GRFDist::computeGeneralizedDistances(trees, metric_spi, ABSOLUTE).number_of_unique_trees, 0);
  EXPECT_EQ(GRFDist::computeGeneralizedDistances(trees, metric_mci, ABSOLUTE).number_of_unique_trees, 0);
  EXPECT_EQ(GRFDist::computeGeneralizedDistances(trees, metric_msi, RELATIVE).number_of_unique_trees, 0);
  EXPECT_EQ(GRFDist::computeGeneralizedDistances(trees, metric_spi, RELATIVE).number_of_unique_trees, 0);
  EXPECT_EQ(GRFDist::computeGeneralizedDistances(trees, metric_mci, RELATIVE).number_of_unique_trees, 0);
}

TEST_F(DistancesTest, ExampleFromSlideshow) {
  PllSplit::setTipCount(6);
  PllTree tree = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  std::vector<PllTree> trees = {tree, tree};
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}
TEST_F(DistancesTest, ComparisionTree0_2taxa24) {
  std::vector<PllTree> trees = load_24taxa();
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}
TEST_F(DistancesTest, example_24_msi) {
  std::vector<PllTree> trees = load_24taxa();
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}

TEST_F(DistancesTest, example_from_slideshow) {
  PllSplit::setTipCount(6);
  PllTree tree1 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  PllTree tree2 = TreeReader::readTreeFile(current_data_dir + "example_from_slideshow")[0];
  std::vector<PllTree> trees = {tree1, tree2};
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_msi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_mci, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
  EXPECT_NEAR(GRFDist::computeGeneralizedDistances(
    trees, metric_spi, ABSOLUTE).pairwise_distance_mtx[0][0], 0, epsilon);
}

TEST_F(DistancesTest, test_instances) {
  std::vector<std::string> instance_names =
        {"24", "125", "141", "143",
        "148", "150", "218", "350",
        "404", "500"};
        //"628", "714", "885", "994",
        //"1288", "1481", "1512", "1604",
        //"1908", "2000"
  for (std::string name : instance_names){
    std::cout << "heads/" << name << " MSI" << std::endl;
    execute_test_generalized("heads/" + name, metric_msi, SIMILARITY);
    std::cout << "heads/" << name << " SPI" << std::endl;
    execute_test_generalized("heads/" + name, metric_spi, SIMILARITY);
    std::cout << "heads/" << name << " MCI" << std::endl;
    execute_test_generalized("heads/" + name, metric_mci, SIMILARITY);
    std::cout << "heads/" << name << " RF absolute" << std::endl;
    execute_test_rf("heads/" + name, ABSOLUTE);
    std::cout << "heads/" << name << " RF relative" << std::endl;
    execute_test_rf("heads/" + name, RELATIVE);
  }
}
