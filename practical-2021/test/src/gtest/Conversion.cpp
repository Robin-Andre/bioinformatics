#include "gtest/gtest.h"
#include "../../../src/io/RFDataReader.hpp"
#include "../../../src/io/RFDataWriter.hpp"
#include "../../../src/io/Converter.hpp"
#include <filesystem>
#include <dirent.h>

class Conversion : public testing::Test {
protected:
  std::string raxml_path = "../test/res/reference_results/";
  std::string target_path = "../test/res/references_json/";
  std::string matrix_path = "../test/res/R_results/";
  //std::vector<std::string> test_setups = {"heads" , "full"};
  std::vector<std::string> test_setups = {"full"};
  std::vector<std::string> test_sizes = {"24"};
  /*"125", "141", "143", "148",
  "150", "218", "350", "354", "404",
  "500", "628", "714", "885", "994",
  "1288", "1481", "1512", "1604",
   "1908", "2000", "2308", "2554"};*/

  std::vector<std::string> metrics = {"MSI" , "SPI", "MCI"};
  std::vector<std::string> modes = {"SIMILARITY" , "RELATIVE", "ABSOLUTE"};
};

/*TEST_F(Conversion, convert_raxml) {
  std::string metric_name = "RF";
  //RAXMLReader::read needs to be adapted to switch from relative to absolute!
  std::string mode_name = "ABSOLUTE";
  for(std::string setup : test_setups){
    for(std::string size : test_sizes){
      std::string in_path = raxml_path + setup + "/" + size + "/";
      std::cout << in_path << std::endl;
      DIR* dir = opendir(in_path.c_str());
      if (dir) {
        Converter::raxmlToJSON(in_path, target_path + 
            metric_name + "/" + mode_name + "/" + setup + "/" + size, mode_name, std::stoi(size));
      }
    }
  }

}

TEST_F(Conversion, convert_matrix) {
  for(std::string mode : modes){
    for(std::string metric : metrics){
      for(std::string setup : test_setups){
        for(std::string size : test_sizes){
          std::string in_path = matrix_path + metric + "/" + mode + "/" + setup + "/" + size;
          std::cout << in_path << std::endl;
          if (FILE *file = fopen(in_path.c_str(), "r")) {
            Converter::matrixToJSON(in_path, target_path + metric + 
                "/" + mode + "/" + setup + "/" + size, metric, mode);
          }
        }
      }
    }
  }
}*/
