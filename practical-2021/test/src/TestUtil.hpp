#pragma once
#include<vector>
#include "../../../src/datastructures/PllSplits.hpp"
class TestUtil {

public:
  //I am not happy with this construction at all, why can't we have this simply in PLLsplit for starters
  static PllSplit createSplit(const std::vector<size_t>& part1) {
    //if (part1[0] != 0) throw "In every split, 0 must be in Partition 1, hence it must hold that part1[0]==0"; 
    auto split_bits = (pll_split_t)calloc(PllSplit::getSplitLen(), sizeof(pll_split_base_t));
    setBits(split_bits, part1);
    return PllSplit(split_bits);
  }

  static void setBits(pll_split_t split_bits, std::vector<size_t> part1) {
    size_t major_idx;
    size_t minor_idx;
    for(size_t tip : part1){
      major_idx = tip / (sizeof(pll_split_base_t) * 8);
      minor_idx = tip % (sizeof(pll_split_base_t) * 8);
      split_bits[major_idx] |= (1 << minor_idx);
    }
  }

  static PllSplitList createSplitList(std::vector<std::vector<size_t>> part1s){
    std::vector<PllSplit> splits;
    size_t split_len = PllSplit::getSplitLen();
    pll_split_t split_pointer = (pll_split_t) calloc(part1s.size()* split_len, sizeof(pll_split_base_t));
    for (size_t i=0; i<part1s.size(); ++i) {
      setBits(split_pointer + i* split_len, part1s[i]);
      splits.emplace_back(PllSplit(split_pointer + i * split_len));
    }
    std::sort(splits.begin(), splits.end());
    return PllSplitList(splits);
  }
};
