#include "PhylogeneticMathUtils.hpp"
#include "io/RFDataWriter.hpp"
#include "Distances.hpp"
#include "Metric.hpp"
#include "io/TreeReader.hpp"

int main(int argc, char* argv[]) {
    io::Config config = io::parseCommandLineOptions(argc, argv);
    std::cout << "Input File Path: " << config.input_file_path  << std::endl;

    io::IOData result;
    if(config.metric == nullptr) {
      std::cout << "Metric Used: RF\n";
      result = Distances::computeRFDistances(
               TreeReader::readTreeFile(config.input_file_path), ABSOLUTE);
    }
    else {
      std::cout << "Metric Used: " << config.metric->name() <<"\n";
      result = Distances::computeGeneralizedDistances(
               TreeReader::readTreeFile(config.input_file_path), *(config.metric), ABSOLUTE);
    }
    std::cout << result.toString();


    if(config.output_file_path.size() <= 1) {
      std::cout << "The Output was not properly specified. no output file will be written\n";
      config.remove();
      exit(0);
    }


    RAXMLWriter::write(config.output_file_path, result);
    config.remove();
}
