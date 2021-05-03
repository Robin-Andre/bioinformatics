#include "PllSplits.hpp"
#include "PllTree.hpp"

/*  This is an example function. It is _slow_. You should replace it */
size_t PllSplit::popcount() {
  size_t popcount = 0;
  for (size_t index = 0; index < _amount_of_register * splitBitWidth(); ++index) {
    if (bitExtract(index) == 1) { popcount += 1; }
  }
  return popcount;
}

uint32_t PllSplit::bitExtract(size_t bit_index) {
  pll_split_base_t split_part = _split[computeMajorIndex(bit_index)];
  return (split_part & (1u << computeMinorIndex(bit_index))) >>
         computeMinorIndex(bit_index);
}

int PllSplit::compareTo(PllSplit other) const {
  for (size_t i = 0; i < _amount_of_register; ++i) {
    if (_split[i] != other()[i]){
      return (int) (_split[i] > other()[i]?1:-1);
    }
  }
  return 0;
}

PllSplitList::PllSplitList(const PllTree &tree) {
  size_t split_len = (tree.getTipCount() / computSplitBaseSize());
  if ((tree.getTipCount()  % computSplitBaseSize()) > 0) { split_len += 1; }
  auto tmp_splits = pllmod_utree_split_create(
      tree.tree()->vroot, tree.tree()->tip_count, nullptr);

  for (size_t i = 0; i < tree.tree()->tip_count - 3; ++i) {
    _splits.emplace_back(PllSplit(tmp_splits[i], split_len));
  }
  free(tmp_splits);
}

PllSplitList::PllSplitList(const std::vector<PllSplit> &splits) {
  if(splits.size() > 0){
    size_t split_len = splits[0].getAmountOfRegister();
    for (size_t i = 0; i < splits.size(); ++i) {
      auto tmp_splits = (pll_split_t)calloc(split_len, sizeof(pll_split_base_t));
      memcpy(tmp_splits, splits[i](), split_len * sizeof(pll_split_base_t));
      _splits.emplace_back(PllSplit(tmp_splits, split_len));
    }
  }
}



PllSplitList::~PllSplitList() {
  if (!_splits.empty()) { free(_splits[0]()); }
}



PllSplitList PllSplitList::symmetricDifference(const PllSplitList& other) const {
  std::vector<PllSplit> different_splits;
  if (_splits.size() == 0) {
    return PllSplitList(other);
  }
  size_t other_split_count = other.getSplitCount();
  if(other_split_count == 0){
    return PllSplitList(_splits);
  }
  size_t i = 0;
  size_t j= 0;
  while (i < _splits.size() && j < other_split_count){
    int cmp = _splits[i].compareTo(other[j]);
    if (cmp == -1) { //this < other
      different_splits.push_back(_splits[i]);
      ++i;
    } else if (cmp == 1){ //this > other
      different_splits.push_back(other[j]);
      ++j;
    } else { //cmp == 0, this = other
      ++i;
      ++j;
    }
  }
  while (i < _splits.size()) {
    different_splits.push_back(_splits[i]);
    ++i;
  }
  while (j < other_split_count) {
    different_splits.push_back(other[j]);
    ++j;
  }
  return PllSplitList(different_splits);

}
