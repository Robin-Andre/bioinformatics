#pragma once
extern "C" {
#include "libpll/pll.h"
#include "libpll/pll_tree.h"
}
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include <immintrin.h>
#include <float.h>
#include <bitset>
#include <algorithm>
#include <sstream>
#include <string>
#include "../enums.hpp"



class PllTree;
class phylomath;

/*
 * A convenience class for the purposes of doing math on the individual splits.
 * The pll_split_t is a non-owning pointer, so this class does not have a
 * destructor. The functions which are probably needed:
 *
 *  - popcount
 *  - bitextract
 *  - intersect (and)
 *  - union (or)
 *  - set minus (xor)
 *
 * I have already written an example bit extract and popcount, but they are far
 * from optimal. Remember that the underlying type is a *pointer* to an array of
 * pll_split_base_t (which is an unsigned int), because the number of taxa in
 * the could be more than 32 (or even 64). This means that you will need to
 * iterate over the array to compute the correct value for popcount etc.
 */
class PllSplit {
public:
  explicit PllSplit(pll_split_t s);
  PllSplit() {
    PllSplit(static_cast<pll_split_t> (calloc(PllSplit::split_len, sizeof(pll_split_base_t))));
  }
  pll_split_t operator()() const { return _split; }
  pll_split_t operator !() const { return _inverse; }
  friend bool operator == (const PllSplit& p1, const PllSplit& p2);
  friend bool operator < (const PllSplit& p1, const PllSplit& p2);

  size_t   popcount() const;
  uint32_t bitExtract(size_t bit_index) const;
  size_t partitionSizeOf (Partition block) const {
    //assert(splitValid());
    return (block == Block_A) ? size_block_A : size_block_B;
  }

  double h() const {return h_value;}
  double entropy() const {return entropy_value;}
  size_t intersectionSize(const PllSplit& other, Partition partition_this, Partition partition_other) const;
  static size_t intersectionSize(const pll_split_t& first, const pll_split_t& second) {
  size_t count = 0; 
  for(unsigned i = 0; i < PllSplit::getSplitLen(); ++i) {
    count += __builtin_popcount(first[i] & second[i]);
  }
  return count;
}


  std::string toString() const;

  static void setTipCount(size_t val) {
    PllSplit::tip_count = val;

    size_t split_base_size = PllSplit::computSplitBaseSize();
    PllSplit::split_len = (PllSplit::tip_count / split_base_size);
    if (tip_count % split_base_size > 0) { PllSplit::split_len += 1; }
    assert(PllSplit::split_len * split_base_size >= tip_count);

    pll_split_base_t bit_mask = 0;
    size_t offset = val - ((PllSplit::split_len - 1) * split_base_size);
    for(size_t i = 0; i < offset; ++i){
      bit_mask |= (1 << i);
    }
    PllSplit::bitmask_for_unused_bits = bit_mask;


  }
  static size_t getTipCount() {
    return PllSplit::tip_count;
  }

  static size_t getSplitLen() {
    return PllSplit::split_len;
  }

private:
  constexpr size_t splitBitWidth() const {
    return sizeof(pll_split_base_t) * 8;
  }

  constexpr size_t computeMajorIndex(size_t index) const {
    return index / splitBitWidth();
  }

  constexpr size_t computeMinorIndex(size_t index) const {
    return index % splitBitWidth();
  }

  /* Computes the number of bits per split base */
  static inline size_t computSplitBaseSize() {
    return sizeof(pll_split_base_t) * 8;
  }

  //bool splitValid() const;
  //size_t basePopcount(pll_split_base_t count) const;

  pll_split_t _split = nullptr;
  pll_split_t _inverse = nullptr;
  size_t size_block_A;
  size_t size_block_B;
  double h_value;
  double entropy_value;


  static size_t tip_count;
  static size_t split_len;
  static pll_split_base_t bitmask_for_unused_bits;

};
class PllSplitList {
public:
  explicit PllSplitList(const PllTree &tree);
  explicit PllSplitList(const std::vector<size_t> &splits);
  //We can now actually use the default constructor since the map is constantly pushing. 
  //TODO buid nondefault constructor which preallocates elements
  explicit PllSplitList() {

  }
  /* Rule of 5 constructors/destructors */
  ~PllSplitList();
  PllSplitList(const PllSplitList &other) : PllSplitList(other._split_offsets) {
    this->maximum_entropy = other.getMaximumEntropy();
    this->maximum_information_content = other.getMaximumInformationContent();
  }
  PllSplitList(PllSplitList &&other) :
      _split_offsets(std::exchange(other._split_offsets, {})) {
        this->maximum_entropy = other.getMaximumEntropy();
        this->maximum_information_content = other.getMaximumInformationContent();
      }
  PllSplitList &operator=(const PllSplitList &other) {
    return *this = PllSplitList(other);
  }
  PllSplitList &operator=(PllSplitList &&other) {
    std::swap(_split_offsets, other._split_offsets);
    return *this;
  }

  friend bool operator == (const PllSplitList& p1, const PllSplitList& p2);
  size_t operator[](size_t index) const { return _split_offsets[index]; }
  //size_t pos(size_t index) const {return _split_offsets[index];}

  const std::vector<size_t>& getSplits() const {return _split_offsets;}
  size_t getSplitCount() const {return _split_offsets.size();}

  double getMaximumEntropy() const {return maximum_entropy;}
  double getMaximumInformationContent() const {return maximum_information_content;}

  std::string toString() const;
  void push(const PllSplit& split, size_t offset);

private:
  /* This vector now represents the PllSplitList, each element is now an array index of
  the big PllMap. 
  */
  std::vector<size_t> _split_offsets; 
  double maximum_entropy = 0;
  double maximum_information_content = 0;
};
