#pragma once
#include <string>
#include <getopt.h>
#include <iostream>
namespace io {
struct Config {
  std::string input_file_path;
  std::string output_file_path;
};
static Config parseCommandLineOptions(int argc, char* argv[]){ 
    int opt;
    Config config; 
    while((opt = getopt(argc, argv, "i:o:" ))!= -1) {
    std::cout << opt << " " << optarg << "\n";
      switch(opt) {
          //std::cout << opt << " " << optarg << "\n";
          case 'i':
            printf("parameter 'i' specified with the value %s\n", optarg);
            config.input_file_path = optarg;
            break;
          case 'o':
            printf("parameter 'o' specified with the value %s\n", optarg);
            config.output_file_path = optarg;
            break;
          case '?':
            std::cout << "Unknown Para\n";
            exit(1);

      }   
    }
    return config;


}
} // namespace io