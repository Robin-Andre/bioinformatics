#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllTree.hpp"
#include "RFData.hpp"
#include "io/TreeReader.hpp"
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <numeric>
#include <bitset>
#include <algorithm>
#include <cfloat>


class RFDistance {
public:
  RFData computeRF(const std::vector<PllTree>& trees);
private:
  //determines position in linear array, make sure that i < j
size_t arrayPos(size_t i, size_t j, size_t tree_count)  {
    assert(i < tree_count && j < tree_count);
    assert(i < j);
    size_t offset =(i*(2*tree_count-i-1))/2;
    return offset + (j - i - 1);
  }
};
