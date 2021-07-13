#include "PhylogeneticMathUtils.hpp"
#include "io/RFDataWriter.hpp"
#include "GeneralizedRFDistance.hpp"
#include "Metric.hpp"
#include "io/TreeReader.hpp"


int main(int argc, char* argv[]) {
    
    auto time_start = std::chrono::high_resolution_clock::now();
    
    io::Config config = io::parseCommandLineOptions(argc, argv);
    std::cout << "Input File Path: " << config.input_file_path  << std::endl;
    std::cout << "Metric Used: " << config.metric->name() <<"\n";
    io::IOData result = GeneralizedRFDistance::computeGeneralizedDistances(TreeReader::readTreeFile(config.input_file_path), *(config.metric), ABSOLUTE);
    
    
    auto time_end = std::chrono::high_resolution_clock::now();
    std::cout << "Time needed (ms): " << std::to_string((std::chrono::duration<double, std::milli>(time_end - time_start)).count()) << "\n";


    std::string result_text = config.input_file_path + " " + config.metric->name() + " " 
                + std::to_string((std::chrono::duration<double, std::milli>(time_end - time_start)).count()) + " "+ GIT_COMMIT_HASH;
    io::write_benchmark_timing(result_text);

    std::cout << result.toString();

    
    if(config.output_file_path.size() <= 1) {
      std::cout << "The Output was not properly specified. no output file will be written\n";
      config.remove();
      exit(0);
    }


    RAXMLWriter::write(config.output_file_path, result);
    config.remove();
}
