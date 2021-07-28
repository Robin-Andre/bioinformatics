# Pairwise RF-Distances
This tool calculates the pairwise distance of a set of phylogenetic trees based on the metrics found in the [TreeDist R package](https://github.com/ms609/TreeDist/tree/2.0.3) as well as the standard [Robinson-Foulds metric](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)
The implemented metrics are Mutual Cluster Information (MSI), Shared Phylogenetic Information (SPI) and Matching Split Information (MSI).
## Mutual Cluster Information
## Shared Phylogenetic Information
## Matching Split Information
## Installation 
### Requirements
The following software is required to run $our_software_name$
- A c++17 ready compiler such as `g++ (version)` or `clang (version)`
- [Google OR-Tools](https://developers.google.com/optimization/install/cpp) 
- cmake
Install using `make full && cd build && make` 

To build without tests run `make && cd build && make`

The binary file $our_binary_name$ will be located in the folder `bin`

## Running $our_software_name$
 To run $our_software_name$ execute `./$binary_name$ `
 
### Command Line Parameters
- (mandatory) -i *path_to_file* specifies a path to a file with phylogenetic trees in the [Newick format](https://en.wikipedia.org/wiki/Newick_format)
- (optional) -o *path_to_file* specifies an output path. Two files will be written an output and an info file. Not specifying will print the diagonal matrix on the console
- (mandatory) -m (MSI/SPI/MCI/RF) specifies the metric for evaluation

## Experimental Results
![Alt text](https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMCI_10_trees.png "MCI")
![Alt text](https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMCI_100_trees.png "MCI")
![Alt text](https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonSPI_10_trees.png "SPI")
![Alt text](https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonSPI_100_trees.png "SPI")
![Alt text](https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMSI_10_trees.png "MSI")
![Alt text](https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMSI_100_trees.png "MSI")
