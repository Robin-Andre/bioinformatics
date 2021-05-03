#pragma once
#include <string>
#include <unistd.h>
#include <iostream>
namespace io {
struct Config {
  std::string input_file_path;
};
static Config parseCommandLineOptions(int argc, char* argv[]){ 
    int opt;
    Config config; 
    while(opt = getopt(argc, argv, "i:" )!= -1) {
    std::cout << opt << " " << optarg << "\n";
      switch(opt) {
          std::cout << opt << " " << optarg << "\n";
          case 1:
            printf("parameter 'i' specified with the value %s\n", optarg);
            config.input_file_path = optarg;
            break;
          case ':':
            std::cout << "Missing Para\n";
            break;
          case '?':
            std::cout << "Unknown Para\n";
            break;

      }   
    }
    return config;


}
} // namespace io