#include "io/CommandLineOptions.hpp"
#include "RFDistance.hpp"
#include "io/FileReader.hpp"
int main(int argc, char* argv[]) {
    io::Config config = io::parseCommandLineOptions(argc, argv);
    std::cout << config.input_file_path << "Thefile\n";
    RFDistance distance_calculator;
     
    RFData result = distance_calculator.computeRF(config.input_file_path);
    io::writeOutput(result, config);
}