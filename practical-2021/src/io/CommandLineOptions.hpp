#pragma once
#include <string>
#include <getopt.h>
#include <iostream>
#include "../Metric.hpp"
#include <gmp.h>
namespace io {
struct Config {
  std::string input_file_path;
  std::string output_file_path;
  Metric* metric;
  //As long as there is a raw pointer this remove is necessary when config runs out of scope
  void remove() {
    delete metric;
  }
};

static Metric* metricFromString(const std::string& s) {
  if(s == "RF") {
    return new RFMetric();
  }
  else if(s == "SPI") {
    return new SPIMetric();
  }
  else if(s == "MCI") {
    return new MCIMetric();
  }
  else if(s == "MSI") {
    return new MSIMetric();
  }
  std::cerr << "Metric misspelled: " << s <<  " exiting...\n";
  exit(1);
}

static Config parseCommandLineOptions(int argc, char* argv[]){ 
    int opt;
    Config config; 
    std::string metric_name = "RF"; //DEFAULT OPTION
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

} // namespace io