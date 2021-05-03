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
#include <bitset>

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
  PllSplit(pll_split_t s, size_t amount_of_register) : _split{s}, _amount_of_register(amount_of_register) {}

  /*~PllSplit() {
    delete [] _split;
  }
  PllSplit(const PllSplit &other) : _amount_of_register(other._amount_of_register){
    _split = &(malloc(sizeof(pll_split_t)));
    memcpy(_split, other._split, sizeof(pll_split_t));
  }

  Array::Array(const Array &other)
    : m_size(other.m_size)
{
    std::cout << "I am copy constructor" << std::endl;

    m_array = new int[m_size];
    std::copy(other.m_array, other.m_array + m_size, m_array);
}
  PllSplit(PllSplit &&other) :
      _split(std::exchange(other._split, {})), split_len(other.split_len) {}*/

  pll_split_t operator()() const { return _split; }

  size_t   popcount();
  uint32_t bitExtract(size_t bit_index);

  int compareTo(PllSplit other) const;


  void printSplit() const {
    std::cout <<_split << ": ";
    for (size_t i = 0; i < _amount_of_register; ++i){
      std::cout << (std::bitset<32>(_split[i]));
    }
  std::cout << std::endl;
  }

  size_t getAmountOfRegister() const {
    return _amount_of_register;}

private:
  constexpr size_t splitBitWidth() const {
    return sizeof(pll_split_base_t) * 8;
  }

  constexpr size_t computeMajorIndex(size_t index) {
    return index / splitBitWidth();
  }

  constexpr size_t computeMinorIndex(size_t index) {
    return index % splitBitWidth();
  }


  pll_split_t _split;
  size_t _amount_of_register;
};

class PllSplitList {
public:
  PllSplitList(const PllTree &tree);
  PllSplitList(const std::vector<PllSplit> &splits);

  /* Rule of 5 constructors/destructors */
  ~PllSplitList();
  PllSplitList(const PllSplitList &other) : PllSplitList(other._splits) {}
  PllSplitList(PllSplitList &&other) :
      _splits(std::exchange(other._splits, {})) {}
  PllSplitList &operator=(const PllSplitList &other) {
    return *this = PllSplitList(other);
  };
  PllSplitList &operator=(PllSplitList &&other) {
    std::swap(_splits, other._splits);
    return *this;
  };

  PllSplit operator[](size_t index) const { return _splits[index]; }

  size_t getSplitCount() const {return _splits.size();}
  PllSplitList symmetricDifference(const PllSplitList& other) const;


  void printSplits() const {
    std::cout << "-------------------------"<< std::endl;
    for(size_t i = 0; i < _splits.size(); ++i){
      _splits[i].printSplit();
    }
    std::cout << "-------------------------"<< std::endl;
  }



private:
  /* Computes the number of bits per split base */
  constexpr size_t computSplitBaseSize() const {
    return sizeof(pll_split_base_t) * 8;
  }



  size_t computeSplitArraySize() const {
    return _splits.size() == 0 ? 0 : _splits[0].getAmountOfRegister() * _splits.size();
  }
  std::vector<PllSplit> _splits;
};
