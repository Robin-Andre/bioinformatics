#pragma once
#include <string>
#include <getopt.h>
#include <iostream>
#include "../Metric.hpp"

//TODO: @Robin: Add doc if needed
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif
namespace io {
  /* This struct contains the relevant configurations neccessary to run the algorithm
   *
   */
struct Config {
  std::string input_file_path;
  std::string output_file_path;
  GeneralizedMetric* metric;
  //As long as there is a raw pointer this remove is necessary when config runs out of scope
  void remove() {
    delete metric;
  }
};

static GeneralizedMetric* metricFromString(const std::string& s) {
  if(s == "SPI") {
    return new SPIMetric();
  }
  else if(s == "MCI") {
    return new MCIMetric();
  }
  else if(s == "MSI") {
    return new MSIMetric();
  }
  else if(s == "RF") {
    return nullptr;
  }
  std::cerr << "Cannot find Metric [ " << s <<  " ] Metrics are specified via -m (MSI/SPI/MCI) exiting...\n";
  exit(1);
}

static Config parseCommandLineOptions(int argc, char* argv[]){
    int opt;
    Config config;
    std::string metric_name = ""; //DEFAULT OPTION IS INVALID
    while((opt = getopt(argc, argv, "i:o:m:" )) != -1) {
      switch(opt) {
          case 'i':
            config.input_file_path = optarg;
            break;
          case 'o':
            config.output_file_path = optarg;
            break;
          case 'm':
            metric_name = optarg;
            break;
          case '?':
            std::cout << "Unknown Parameter exiting ...\n";
            exit(1);

      }
    }
    config.metric = io::metricFromString(metric_name);
    return config;
}
static void write_benchmark_timing(std::string& input) {
  std::ofstream out_file(std::string("../benchmark/") + GIT_COMMIT_HASH + ".csv", std::ios::app);
  if(out_file.is_open()) {
    out_file << input << "\n";
  }
  else {
    std::cerr << "Benchmark Timing File cannot be opened\n";
  }
  out_file.close();
}
} // namespace io
