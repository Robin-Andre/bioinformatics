#include "PllSplits.hpp"
#include "PllTree.hpp"

size_t PllSplit::split_len = 0;

/*  This is an example function. It is _slow_. You should replace it */
size_t PllSplit::popcount() {
  size_t popcount = 0;
  for (size_t index = 0; index < PllSplit::split_len * splitBitWidth(); ++index) {
    if (bitExtract(index) == 1) { popcount += 1; }
  }
  return popcount;
}
/*@Luise This is the stuff I am really not proud of :( The two operators are needed for
sorting and they are far from correct. Right now they only check the first register instead of all
of them and I have no clue how we are gonna pass the information of how many registers are actually
required for the corresponding amount of taxa. The information is there
(for example in the PllSplitlist.computeSplitLen())
but I lack a good idea to integrate it. What is worse that the < operator might even be correct for most cases
the == operator is essentially always wrong for taxa > 32. In the end it might even turn out that the
Split infer from plllib is already sorted and that all of this is void.

Also I have no idea how the taxa are stored. Is it MSB or LSB? Take a 33 taxa tree for example.
Would the first taxa X be either located :
          Register[0]                  |            Register[1]
00000000000000000000000000000000       | 0000000000000000000000000000000X <- here
0000000000000000000000000000000X <-here| 00000000000000000000000000000000
Not that it matters for sorting at all.  just that I am clueless

Also if the two operators are working properly then so should the (currently inefficient) algorithm
  */
bool operator == (const PllSplit& p1, const PllSplit& p2) {
  for(unsigned i = 0; i < PllSplit::split_len; ++i) {
    if(p1._split[i] != p2._split[i]) {
      return false;
    }
  }
  return true;
  //return p1._split[0] == p2._split[0] && p1._split[1] == p2._split[1]&& p1._split[2]
  //== p2._split[2] && p1._split[3] == p2._split[3];
  //The way to fix this would be r[0] == s[0] && r[1]==s[1] && .. &&r[n]==s[n] but the ominous number n is missing
}
bool operator < (const PllSplit&p1, const PllSplit& p2) {
  for(unsigned i = 0; i < PllSplit::split_len; ++i) {
    if(p1._split[i] != p2._split[i]) {
      return (p1._split[i] < p2._split[i]);
    }
  }
  return false;
  //return !(p1._split[0] >= p2._split[0] && p1._split[1] >= p2._split[1]
  //&& p1._split[2] >= p2._split[2] && p1._split[3] >= p2._split[3]); /
  // Similar to above, some cool way to fix it if the info of registers is known
}

uint32_t PllSplit::bitExtract(size_t bit_index) const {
  assert(bit_index < PllSplit::split_len * computSplitBaseSize());
  pll_split_base_t split_part = _split[computeMajorIndex(bit_index)];
  return (split_part & (1u << computeMinorIndex(bit_index))) >> computeMinorIndex(bit_index);
}

PllSplitList::PllSplitList(const PllTree &tree) {
  assert(PllSplit::getSplitLen() == PllSplit::computeSplitLen(tree.getTipCount()));
  pll_split_t* tmp_splits = pllmod_utree_split_create(tree.tree()->vroot, tree.getTipCount(), nullptr);
  //_splits.reserve(tree.tree()->tip_count - 3);
  for (size_t i = 0; i < tree.getTipCount() - 3; ++i) {
    _splits.emplace_back(PllSplit(tmp_splits[i]));
  }
  free(tmp_splits);
}

PllSplitList::PllSplitList(const std::vector<PllSplit> &splits) {
  if(splits.size() > 0){
    pll_split_t split_pointer = (pll_split_t) calloc(splits.size()* PllSplit::getSplitLen(), sizeof(pll_split_base_t));
    for (size_t i=0; i<splits.size(); ++i) {
      memcpy(split_pointer + i*PllSplit::getSplitLen(), splits[i](), PllSplit::getSplitLen() * sizeof(pll_split_base_t));
      _splits.emplace_back(PllSplit(split_pointer + i*PllSplit::getSplitLen()));
    }
  }
}



PllSplitList PllSplitList::symmetricDifference(const PllSplitList& other) const {
  std::vector<PllSplit> different_splits;
  if (_splits.size() == 0) return PllSplitList(other);
  size_t other_split_count = other.getSplitCount();
  if(other_split_count == 0) return PllSplitList(_splits);
  size_t i = 0;
  size_t j= 0;
  while (i < _splits.size() && j < other_split_count){
    if (_splits[i] == other[j]) {
      ++i;
      ++j;
    } else if (_splits[i] < other[j]) {
      different_splits.push_back(_splits[i]);
      ++i;
    } else {
      different_splits.push_back(other[j]);
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

size_t PllSplitList::rfDistance(const PllSplitList& other) const {
  size_t other_split_count = other.getSplitCount();
  if (_splits.size() == 0) return other_split_count;
  if(other_split_count == 0) return _splits.size();
  size_t i = 0;
  size_t j = 0;
  size_t distance = 0;
  while (i < _splits.size() && j < other_split_count){
    if (_splits[i] == other[j]) {
      ++i;
      ++j;
    } else if (_splits[i] < other[j]) {
      ++distance;
      ++i;
    } else {
      ++distance;
      ++j;
    }
  }
  distance += (_splits.size() - i);
  distance += (other_split_count - j);
  return distance;
}

PllSplitList::~PllSplitList() {
  if (!_splits.empty()) { free(_splits[0]()); }
}

bool operator == (const PllSplitList& p1, const PllSplitList& p2) {
  std::vector<PllSplit> splits1 = p1.getSplits();
  std::vector<PllSplit> splits2 = p2.getSplits();
  if(splits1.size() != splits2.size()) return false;
  for(size_t i = 0; i < splits1.size(); ++i){
    if (!(splits1[i] == splits2[i])) return false;
  }
}
