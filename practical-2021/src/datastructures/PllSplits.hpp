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

class PllTree;
class phylomath;

/*
 * Represents a split, which corresponds to a branch in a phylogenetic tree.
 * If the respective branch is removed from the tree,
 * two disconnected subtrees are produced.
 * Depending on the subtree, which contains a certain taxon,
 * the taxa are assigend to the partitions of the split.
 * According to that, bits are set in the bitvector representation
 * which is encapsulated by this class.
 */

typedef bool partition_t;
class PllSplit {
public:
  /**
   * TODO: Robin, please help here!
   */
  explicit PllSplit(pll_split_t s);
  PllSplit() {
    size_block_A = 0;
    size_block_B = 0;
    h_value = 0.0;
    entropy_value = 0.0;
    PllSplit(static_cast<pll_split_t> (calloc(PllSplit::split_len, sizeof(pll_split_base_t))));
  }
  pll_split_t operator()() const { return _split; }
  friend bool operator == (const PllSplit& p1, const PllSplit& p2);
  friend bool operator < (const PllSplit& p1, const PllSplit& p2);

  /**
   * Gives the size of the partition of the PllSplit represented by 1 or 0 resp.
   *
   * @param block: partition (1 or 0)
   * @return: size of the partition
   */
  size_t partitionSizeOf (partition_t block) const {
    return block ? size_block_A : size_block_B;
  }

  /**
   * @return: The information content of the PllSplit (see Phylomath for details)
   */
  double h() const {return h_value;}

  /**
   * @return: The entropy of the PllSplit (see Phylomath for details)
   */
  double entropy() const {return entropy_value;}

  /**
   * Determines the size of the intersection of the 1-partition of this PllSplit
   * and the 1-partition of the provided PllSplit
   *
   * @param other: split to intersect this PllSplit with
   * @return: size of the intersection
   */
  size_t intersectionSize(const PllSplit& other) const;

  /**
   * @return: A string representation for the PllSplit
   */
  std::string toString() const;


  /**
   * Globally sets the tip count for all PllSplits and updates other values
   * depending on the tip count
   *
   * @param val: The new tip count
   */
  static void setTipCount(size_t val) {
    PllSplit::tip_count = val;

    size_t split_base_size = sizeof(pll_split_base_t) * 8;
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

  /**
   * @return: Global Tip Count
   */
  static size_t getTipCount() {
    return PllSplit::tip_count;
  }

  /**
   * @return: The number of registers required to store the PllSplit
   */
  static size_t getSplitLen() {
    return PllSplit::split_len;
  }

private:
  size_t  popcount() const;


  pll_split_t _split = nullptr;
  size_t size_block_A;
  size_t size_block_B;
  double h_value;
  double entropy_value;


  static size_t tip_count;
  static size_t split_len;
  static pll_split_base_t bitmask_for_unused_bits;

};


/*
 * Encapsulates all splits corresponding to the branches in
 * a phylogenetic tree.
 */

class PllSplitList {
public:
  /**
   * Creates a split list representing the given tree
   *
   * @param tree: Tree for which a split list is to be created
   */
  explicit PllSplitList(const PllTree &tree);

  /**
   * Creates a split list from the provided vector of splits
   *
   * @param splits: splits to encapsulate in a SplitList
   */
  explicit PllSplitList(const std::vector<size_t> &splits);


  /**
   * TODO: Robin please help here
   */
  //TODO buid nondefault constructor which preallocates elements
  explicit PllSplitList() {

  }
  /* Rule of 5 constructors/destructors */
  ~PllSplitList(){}
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

  /**
   * Creates a PllSplitList from the provided vector of PllSplits
   *
   * @return The splits encapsulated in the Split List Object as vector
   */
  const std::vector<size_t>& getSplits() const {return _split_offsets;}

  /**
   * @return The number of PllSplits in the PllSplitList
   */
  size_t getSplitCount() const {return _split_offsets.size();}

  /**
   * @return The maximum entropy, i.e. the sum of the entropy of all PllSplits in
   * the PllSplitList
   */
  double getMaximumEntropy() const {return maximum_entropy;}

  /**
   * @return The maximum information content , i.e. the sum of the entropy
   * of all splits in the list
   */
  double getMaximumInformationContent() const {return maximum_information_content;}

  /**
   * @return: A string representation for the PllSplitList
   */
  std::string toString() const;

//TODO: @Robin: Is this correct?
  /**
   * Inserts a new PllSplit to the PllSplitList
   *
   * @param split: The PllSplit to be inserted
   * @param offset: The position where the PllSplit is to be inserted
   *
   */
  void push(const PllSplit& split, size_t offset);

private:
  /* This vector now represents the PllSplitList, each element is now an array index of
  the big PllMap.
  */
  std::vector<size_t> _split_offsets;
  double maximum_entropy = 0;
  double maximum_information_content = 0;
};
