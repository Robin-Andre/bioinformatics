#include "PllSplits.hpp"
#include "PllTree.hpp"
#include "../enums.hpp"

size_t PllSplit::tip_count = 0;

/*  This is an example function. It is _slow_. You should replace it */
size_t PllSplit::popcount() const{
  assert(splitValid());
  size_t popcount = 0;
  size_t split_len = PllSplit::getSplitLen();
  for(size_t i = 0; i < split_len; ++i){
    popcount+=basePopcount(_split[i]);
  }
  return popcount;
}

bool operator == (const PllSplit& p1, const PllSplit& p2) {
  size_t split_len = PllSplit::getSplitLen();
  for(unsigned i = 0; i < split_len; ++i) {
    if(p1._split[i] != p2._split[i]) {
      return false;
    }
  }
  return true;
}
bool operator < (const PllSplit&p1, const PllSplit& p2) {
  size_t split_len = PllSplit::getSplitLen();
  for(unsigned i = 0; i < split_len; ++i) {
    if(p1._split[i] != p2._split[i]) {
      return (p1._split[i] < p2._split[i]);
    }
  }
  return false;
}
//TODO, we should kill this once we are certain that we don't need it
uint32_t PllSplit::bitExtract(size_t bit_index) const {
  assert(splitValid());
  //assert(bit_index < PllSplit::getTipCount()); //TODO reenable this assertion if the test for validity is deprecated
  pll_split_base_t split_part = _split[computeMajorIndex(bit_index)];
  return (split_part & (1u << computeMinorIndex(bit_index))) >> computeMinorIndex(bit_index);
}
//TODO inline and move to header
size_t PllSplit::partitionSizeOf (Partition block) const {
  assert(splitValid());
  return (block == Block_A) ? this->popcount() : PllSplit::getTipCount() - this->popcount();
}

//TODO think about rewriting this/other to A/B?
size_t PllSplit::intersectionSize(const PllSplit& other, partition_t partition_this, partition_t partition_other) const {
  assert(splitValid());
  assert(other.splitValid());
  size_t split_len = PllSplit::getSplitLen();
  pll_split_base_t this_mask = partition_this ? 0 : ~0; //TODO explanatory text
  pll_split_base_t other_mask = partition_other ? 0 : ~0;
  size_t count = 0;

  for (size_t i = 0; i < split_len - 1; ++i){
    count += basePopcount((_split[i] ^ this_mask) & (other()[i] ^ other_mask));
  }
  count += basePopcount((_split[split_len - 1] ^ this_mask) & (other()[split_len - 1] ^ other_mask)
                        & bitmaskForUnusedBits());
  return count;
}


//TODO @Robin might wanna do a speedtest, or find another implementation with registers
size_t PllSplit::basePopcount(pll_split_base_t val) const {
  return std::bitset<32>(val).count();
}

//DONE find out why this is needed / why plllib does fail
/* Robin: It seems that the second condition is failing when the PllSplitList Object is destroyed
   but a copy of the underlying splits i.e. a vector of splits has been made and it worked on.
   One theory is pointer management. (or the lack thereof)
*/
//TODO: Apply bitmask to last register
bool PllSplit::splitValid() const {
  //This condition sometimes fails on the splits returned from pll lib, needs to be examined!
  //return (_split != nullptr) && !(_split[0] & ~bitmaskForUnusedBits()) && _split[0] & 1u;
  return (_split != nullptr) &&  _split[0] & 1u;
}
//TODO: Save as field, set together with tipcount
pll_split_base_t PllSplit::bitmaskForUnusedBits() const {
  pll_split_base_t bit_mask = 0;
  size_t offset = PllSplit::getTipCount() - ((PllSplit::getSplitLen() - 1) * computSplitBaseSize());
  //TODO@ Robin, do one shift instead of a loop
  for(size_t i = 0; i < offset; ++i){
    bit_mask |= (1 << i);
  }
  //bit_mask = (pll_split_base_t) (std::pow(2, offset) - 1); //This is hacky but might work
  return bit_mask;
}


//TODO find out if preallocation can be done.
PllSplitList::PllSplitList(const PllTree &tree) {
  assert(PllSplit::getTipCount() == tree.getTipCount());
  pll_split_t* tmp_splits = pllmod_utree_split_create(tree.tree()->vroot, tree.getTipCount(), nullptr);
  //_splits.reserve(tree.tree()->tip_count - 3);
  for (size_t i = 0; i < tree.getTipCount() - 3; ++i) {
    _splits.emplace_back(PllSplit(tmp_splits[i]));
  }
  free(tmp_splits);
}
//TODO find out if preallocation can be done
PllSplitList::PllSplitList(const std::vector<PllSplit> &splits) {
  size_t split_len = PllSplit::getSplitLen();
  if(splits.size() > 0){
    pll_split_t split_pointer = (pll_split_t) calloc(splits.size()* split_len, sizeof(pll_split_base_t));
    for (size_t i=0; i<splits.size(); ++i) {
      memcpy(split_pointer + i* split_len, splits[i](), split_len * sizeof(pll_split_base_t));
      _splits.emplace_back(PllSplit(split_pointer + i*split_len));
    }
  }
}

PllSplitList::~PllSplitList() {
  if (!_splits.empty()) { free(_splits[0]()); }
}
bool operator == (const PllSplitList& p1, const PllSplitList& p2) {
  const std::vector<PllSplit>& splits1 = p1.getSplits();
  const std::vector<PllSplit>& splits2 = p2.getSplits();
  if(splits1.size() != splits2.size()) return false;
  for(size_t i = 0; i < splits1.size(); ++i){
    if (!(splits1[i] == splits2[i])) return false;
  }
  return true;
}
