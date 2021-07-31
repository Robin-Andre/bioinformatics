# Pairwise RF-Distances
This tool calculates the pairwise distance of a set of phylogenetic trees based on the metrics found in the [TreeDist R package](https://github.com/ms609/TreeDist/tree/2.0.3) as well as the standard [Robinson-Foulds metric](https://www.sciencedirect.com/science/article/abs/pii/0025556481900432?via%3Dihub)
The implemented metrics are Mutual Cluster Information (MSI), Shared Phylogenetic Information (SPI) and Matching Split Information (MSI).

## Implementation Gist
The software calculates the pairwise distances in 3 steps: 

1. Making splits unique: Many distinct trees contain the same split. To reduce calculation time the trees are scanned and mapped to a ordered list of unique splits.

<img align="right" width="60%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/animation1_practical.gif">  


2. The pairwise unique splits are evaluated on the underlying metric (MSI, SPI, MCI) and stored in a global table which is calculated once. This reduces the amount of necessary calculations by up to 40% depending on the instance. Most beneficial are instances with a high amount of very similar trees as the amount of unique splits is relatively small.

3. The matchings between trees is calculated in parallel using OR-tools. Since the split scores of the metric are immutable and already precalculated the parallelization can be done without hassle. The matching step requires >70% runtime even in parallel and provides ample opportunities for future optimization.

## Installation 
### Requirements
The following software is required to run $our_software_name$
- A c++17 ready compiler such as `g++ > 6.0` or `clang > 5.0`
- [Google OR-Tools](https://developers.google.com/optimization/install/cpp) 
- cmake > 3.10  

Install using `make full && cd build && make` 

To build without tests run `make && cd build && make`

The binary file will be located in the folder `bin/`


### Command Line Parameters
- (mandatory) -i *path_to_file* specifies a path to a file with phylogenetic trees in the [Newick format](https://en.wikipedia.org/wiki/Newick_format)
- (optional) -o *path_to_file* specifies an output path. Two files will be written an output and an info file.
- (mandatory) -m (MSI/SPI/MCI/RF) specifies the metric for evaluation

### Example Calls
 To run an example call just copy and paste the following code in the `bin/` folder.
 
 `./rfdist -i ../test/res/data/heads/24 -m MSI` without output files  or 
 
 `./rfdist -i ../test/res/data/heads/24 -m MCI -o ../foo/` with output files

## Code Quality
TODO SOFTWIPE SCORE

## Experimental Results
The Experiments have been performed on Ubuntu 20.04 with a AMD Ryzen 5 2500U Radeon Vega Mobile Gfx @2.0Ghz and L1 128KiB, L2 2MiB, L3 4 MiB, 8GB RAM
The software was compiled via installation guide using g++ 10.3. The TreeDist R Package was installed via the R installer. The dataset can be found [here](https://cme.h-its.org/exelixis/resource/download/software/data.tbz).

The experiments have been run on the first 10/100 trees of the dataset for each of the three new metrics respectively. 

<img align="left" width="47%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMSI_10_trees.png">
<img align="left" width="47%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMSI_100_trees.png">

<img align="left" width="47%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonSPI_10_trees.png">
<img align="left" width="47%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonSPI_100_trees.png">


<img align="left" width="47%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMCI_10_trees.png">
<img align="left" width="47%" src="https://github.com/Robin-Andre/bioinformatics/blob/main/practical-2021/plots_python/output/reference_comparisonMCI_100_trees.png">



