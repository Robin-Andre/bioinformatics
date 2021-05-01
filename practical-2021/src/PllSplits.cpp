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
/*@Luise This is the stuff I am really not proud of :( The two operators are needed for
sorting and they are far from correct. Right now they only check the first register instead of all
of them and I have no clue how we are gonna pass the information of how many registers are actually
required for the corresponding amount of taxa. The information is there (for example in the PllSplitlist.computeSplitLen())
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
  return p1._split[0] == p2._split[0]; //The way to fix this would be r[0] == s[0] && r[1]==s[1] && .. &&r[n]==s[n] but the ominous number n is missing
}
bool operator < (const PllSplit&p1, const PllSplit& p2) {
  return p1._split[0] < p2._split[0]; // Similar to above, some cool way to fix it if the info of registers is known
}

uint32_t PllSplit::bitExtract(size_t bit_index) const {
  pll_split_base_t split_part = _split[computeMajorIndex(bit_index)];
  return (split_part & (1u << computeMinorIndex(bit_index))) >>
         computeMinorIndex(bit_index);
}
PllSplitList::PllSplitList(const PllTree &tree) {
  auto tmp_splits = pllmod_utree_split_create(
      tree.tree()->vroot, tree.tree()->tip_count, nullptr);

  for (size_t i = 0; i < tree.tree()->tip_count - 3; ++i) {
    _splits.emplace_back(tmp_splits[i]);
  }
  free(tmp_splits);
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
