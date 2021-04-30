#include "PllSplits.hpp"
#include "PllTree.hpp"

/*  This is an example function. It is _slow_. You should replace it */
size_t PllSplit::popcount(size_t len) {
  size_t popcount = 0;
  for (size_t index = 0; index < len * splitBitWidth(); ++index) {
    if (bitExtract(index) == 1) { popcount += 1; }
  }
  return popcount;
}

uint32_t PllSplit::bitExtract(size_t bit_index) {
  pll_split_base_t split_part = _split[computeMajorIndex(bit_index)];
  return (split_part & (1u << computeMinorIndex(bit_index))) >>
         computeMinorIndex(bit_index);
}

int PllSplit::compareTo(PllSplit other, size_t split_len) const {
  pll_split_t other_split = other();
  unsigned int i;
  for (i = 0; i < split_len; i++) {
    if (_split[i] != other_split[i]){
      return (int) (_split[i] > other_split[i]?1:-1);
    }
  }
  return 0;
}

PllSplitList::PllSplitList(const PllTree &tree) {
  auto tmp_splits = pllmod_utree_split_create(
      tree.tree()->vroot, tree.tree()->tip_count, nullptr);

  for (size_t i = 0; i < tree.tree()->tip_count - 3; ++i) {
    _splits.emplace_back(tmp_splits[i]);
  }
  free(tmp_splits);
  size_t tip_count = _splits.size() + 3;
  this->split_len = (tip_count / computSplitBaseSize());
  if ((tip_count % computSplitBaseSize()) > 0) { split_len += 1; }
}

PllSplitList::PllSplitList(const std::vector<PllSplit> &splits, size_t len) {
  split_len = len;
  auto tmp_splits = (pll_split_t)calloc(splits.size() * split_len,
                                        sizeof(pll_split_base_t));
  if(splits.size() > 0){
    memcpy(tmp_splits,
           splits[0](),
           splits.size() * split_len * sizeof(pll_split_base_t));
    for (size_t i = 0; i < splits.size() * split_len; ++i) {
      _splits.emplace_back(tmp_splits + (i * split_len));
    }
  }
}

PllSplitList::PllSplitList(const PllSplitList &other) {
  this->split_len = other.split_len;
  auto tmp_splits = (pll_split_t)calloc(other.computeSplitArraySize(),
                                        sizeof(pll_split_base_t));
  if(other._splits.size() > 0){
    memcpy(tmp_splits,
           other._splits[0](),
           other.computeSplitArraySize() * sizeof(pll_split_base_t));

    for (size_t i = 0; i < other.computeSplitArraySize(); ++i) {
      _splits.emplace_back(tmp_splits + (i * other.split_len));
    }
  }

}

PllSplitList::~PllSplitList() {
  if (!_splits.empty()) { free(_splits[0]()); }
}



PllSplitList PllSplitList::symmetricDifference(const PllSplitList& other) const {
  //assert(split_len == other.getSplitLen());
  std::vector<PllSplit> different_splits;
  if (_splits.size() == 0) {
    return PllSplitList(other);
  }
  size_t other_split_count = other.getSplitCount();
  if(other_split_count == 0){
    return PllSplitList(_splits, split_len);
  }
  size_t i = 0;
  size_t j= 0;
  while (i < _splits.size() && j < other_split_count){
    int cmp = _splits[i].compareTo(other[j], split_len);
    if (cmp == -1) { //this < other
      different_splits.push_back(_splits[i]);
      i++;
    } else if (cmp == 1){ //this > other
      different_splits.push_back(other[j]);
      j++;
    } else { //cmp == 0, this = other
      i++;
      j++;
    }
  }
  while (i < _splits.size()) {
    different_splits.push_back(_splits[i]);
    i++;
  }
  while (j < other_split_count) {
    different_splits.push_back(other[j]);
    j++;
  }
  return PllSplitList(different_splits, split_len);

}
