#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include "PllSplits.hpp"
#include "PllTree.hpp"
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
namespace measure{
size_t rf_distance(const PllSplitList& p1, const PllSplitList& p2);
size_t rf_distance(const PllTree& t1, const PllTree& t2);
std::vector<size_t> full_calculation(const std::string& file);
} //namespace measure