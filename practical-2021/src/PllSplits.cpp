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
  pll_split_t other_split = other.getSplit();
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
}

PllSplitList::PllSplitList(const std::vector<PllSplit> &splits) {
  _splits=splits;
}

PllSplitList::PllSplitList(const PllSplitList &other) {
  auto tmp_splits = (pll_split_t)calloc(other.computeSplitArraySize(),
                                        sizeof(pll_split_base_t));

  memcpy(tmp_splits,
         other._splits[0](),
         other.computeSplitArraySize() * sizeof(pll_split_base_t));

  for (size_t i = 0; i < other.computeSplitArraySize(); ++i) {
    _splits.emplace_back(tmp_splits + (i * other.computeSplitLen()));
  }
}

PllSplitList::~PllSplitList() {
  if (!_splits.empty()) { free(_splits[0]()); }
}



PllSplitList* PllSplitList::symmetricDifference(PllSplitList* other) const {
  std::vector<PllSplit> different_splits;
  size_t i = 0;
  size_t j= 0;
  size_t split_len = computeSplitLen();
  size_t other_split_count = other->getSplitCount();
  while (i < _splits.size() && j < other_split_count){
    int cmp = _splits[i].compareTo((*other)[j], split_len);
    if (cmp == -1) { //this < other
      different_splits.push_back(_splits[i]);
      i++;
    } else if (cmp == 1){ //this > other
      different_splits.push_back(_splits[j]);
      j++;
    } else { //cmp == 0, this = other
      i++;
      j++;
    }
  }
  return new PllSplitList(different_splits);

}
