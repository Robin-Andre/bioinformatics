#include "io/CommandLineOptions.hpp"
#include "RFDistance.hpp"
#include "io/RFDataWriter.hpp"
int main(int argc, char* argv[]) {
    io::Config config = io::parseCommandLineOptions(argc, argv);
    RFDistance distance_calculator;
    RFData result = distance_calculator.computeRF(TreeReader::readTreeFile(config.input_file_path));
    if(config.output_file_path.size() <= 1) {
      std::cerr << "The Output was not properly specified...exiting";
      exit(1);
    }
    RAXMLWriter::write(config.output_file_path, result);
}
