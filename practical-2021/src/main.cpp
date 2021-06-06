#include "PhylogeneticMathUtils.hpp"
#include "io/CommandLineOptions.hpp"
#include "io/RFDataWriter.hpp"
#include "GeneralizedRFDistance.hpp"
#include "Metric.hpp"
#include "io/TreeReader.hpp"
#include "io/IOData.hpp"

//TODO @Robin, adapt so that it calls Generalized, or switches between Gen/Standard
int main(int argc, char* argv[]) {
    io::Config config = io::parseCommandLineOptions(argc, argv);
    std::cout << "Input File Path: " << config.input_file_path  << std::endl;
    std::cout << "Metric Used: " << config.metric->name() <<"\n";
    io::IOData result = GeneralizedRFDistance::computeDistances(TreeReader::readTreeFile(config.input_file_path), *(config.metric), false);
    //TODO This stuff should be in the Writer, not in main
    if(config.output_file_path.size() <= 1) {
      std::cerr << "The Output was not properly specified...exiting";
      exit(1);
    }
    RAXMLWriter::write(config.output_file_path, result);
    config.remove();
}
