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
#include <bitset>
#include <algorithm>


class PllTree;

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
  PllSplit(pll_split_t s) : _split{s} {}
  PllSplit() {
    PllSplit((pll_split_t)calloc(PllSplit::getSplitLen(), sizeof(pll_split_base_t)));
  }

  pll_split_t operator()() const { return _split; }
  size_t   popcount();
  uint32_t bitExtract(size_t bit_index) const;
  //pll_split_t getSplit() {return _split;}

  friend bool operator == (const PllSplit& p1, const PllSplit& p2);
  friend bool operator < (const PllSplit& p1, const PllSplit& p2);

  static void setTipCount(size_t val) {
    PllSplit::tip_count = val;
  }


  static size_t getTipCount() {
    return PllSplit::tip_count;
  }

  static size_t getSplitLen() {
    size_t split_len = (PllSplit::tip_count / computSplitBaseSize());
    if (tip_count % computSplitBaseSize() > 0) { split_len += 1; }
    assert(split_len * computSplitBaseSize() >= tip_count);
    return split_len;
  }

  //void intersect(const PllSplit& other, PllSplit result, bool invert_this, bool invert_other) const;
  size_t intersectcount(const PllSplit& other, bool invert_this, bool invert_other) const;
  void invert(PllSplit result) const;



  void printSplit() const {
    std::cout << this << ": "<<_split << ": ";
    size_t split_len = PllSplit::getSplitLen();
    for (size_t i = 0; i < split_len; ++i){
      auto str = std::bitset<32>(_split[i]).to_string();
      std::reverse(str.begin(), str.end());
      std::cout << str << "|";
    }
  std::cout << std::endl;
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
  static size_t computSplitBaseSize() {
    return sizeof(pll_split_base_t) * 8;
  }

  pll_split_base_t bitmaskForUnusedBits() const;
  size_t basePopcount(pll_split_base_t count) const;

  pll_split_t _split = nullptr;
  static size_t tip_count;
};
//bool operator == (const PllSplit & p1, const PllSplit& p2);
class PllSplitList {
public:
  explicit PllSplitList(const PllTree &tree);
  explicit PllSplitList(const std::vector<PllSplit> &splits);

  /* Rule of 5 constructors/destructors */
  ~PllSplitList();
  PllSplitList(const PllSplitList &other) : PllSplitList(other._splits) {}
  PllSplitList(PllSplitList &&other) :
      _splits(std::exchange(other._splits, {})) {}
  PllSplitList &operator=(const PllSplitList &other) {
    return *this = PllSplitList(other);
  }
  PllSplitList &operator=(PllSplitList &&other) {
    std::swap(_splits, other._splits);
    return *this;
  }

  friend bool operator == (const PllSplitList& p1, const PllSplitList& p2);
  PllSplit operator[](size_t index) const { return _splits[index]; }

  std::vector<PllSplit> getSplits() const {return _splits;}
  size_t getSplitCount() const {return _splits.size();}
  PllSplitList symmetricDifference(const PllSplitList& other) const;
  size_t rfDistance(const PllSplitList& other) const;

  /* Reenable for print debugging
  void printSplits() const {
    std::cout << "-------------------------"<< std::endl;
    for(size_t i = 0; i < _splits.size(); ++i){
      _splits[i].printSplit();
    }
    std::cout << "-------------------------"<< std::endl;
  }
  */


private:
  std::vector<PllSplit> _splits;
};
