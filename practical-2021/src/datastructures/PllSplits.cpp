#include "PllSplits.hpp"
#include "PllTree.hpp"

size_t PllSplit::tip_count = 0;

/*  This is an example function. It is _slow_. You should replace it */
size_t PllSplit::popcount() const{
  assert(splitValid());
  size_t popcount = 0;
  size_t split_len = PllSplit::getSplitLen();
  //TODO here we can also move the i = 0 call out of the loop
  for(size_t i = 0; i < split_len; ++i){
    if (i == 0){
      popcount+=basePopcount(_split[i] & bitmaskForUnusedBits());
    } else {
      popcount+=basePopcount(_split[i]);
    }

  }
  /*for (size_t index = 0; index < PllSplit::getTipCount(); ++index) {
    if (bitExtract(index) == 1) { popcount += 1; }
  }*/
  return popcount;
}

bool operator == (const PllSplit& p1, const PllSplit& p2) {
  size_t split_len = PllSplit::getSplitLen();
  //TODO this NEEDS the check for unused Bit mask
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
    //TODO this needs the unused Bit mask IF we use it
    if(p1._split[i] != p2._split[i]) {
      return (p1._split[i] < p2._split[i]);
    }
  }
  return false;
}
//TODO, we should kill this once we are certain that we don't need it
uint32_t PllSplit::bitExtract(size_t bit_index) const {
  assert(splitValid());
  //assert(bit_index < PllSplit::getTipCount());
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
  //TODO Put i == 0 here; and let the loop start from 1 
  for (size_t i = 0; i < split_len; ++i){
    if (i == 0){
      count += basePopcount(((_split[i] ^ this_mask) & (other()[i] ^ other_mask)) & bitmaskForUnusedBits());
    } else {
      count += basePopcount((_split[i] ^ this_mask) & (other()[i] ^ other_mask));
    }
  }
  return count;
}
//TODO @Luise this can be removed but we are going to leave it for SPI
size_t PllSplit::unionSize(const PllSplit& other, partition_t partition_this, partition_t partition_other) const {
  assert(splitValid());
  assert(other.splitValid());
  size_t split_len = PllSplit::getSplitLen();
  pll_split_base_t this_mask = partition_this ? 0 : ~0;
  pll_split_base_t other_mask = partition_other ? 0 : ~0;
  size_t count = 0;
 
  for (size_t i = 0; i < split_len; ++i){
    if (i == 0){ 
      count += basePopcount(((_split[i] ^ this_mask) | (other()[i] ^ other_mask)) & bitmaskForUnusedBits());
    } else {
      count += basePopcount((_split[i] ^ this_mask) | (other()[i] ^ other_mask));
    }
  }
  return count;
}
//TODO @Luise can be removed...... aaaaah maybe we need it
bool PllSplit::containsAsSubset(const PllSplit& other, partition_t partition_this, partition_t partition_other) const {
  assert(splitValid());
  assert(other.splitValid());
  size_t split_len = PllSplit::getSplitLen();
  pll_split_base_t this_mask = partition_this ? 0 : ~0;
  pll_split_base_t other_mask = partition_other ? 0 : ~0;
  size_t count = 0;
  for (size_t i = 0; i < split_len; ++i){
    if (i == 0){
      if ((((_split[i] ^ this_mask) & bitmaskForUnusedBits())| ((other()[i] ^ other_mask) & bitmaskForUnusedBits())) != ((_split[i] ^ this_mask) & bitmaskForUnusedBits())) {
        return false;
      }
    } else {
      if ((_split[i] ^ this_mask) | (other()[i] ^ other_mask) != (_split[i] ^ this_mask)) {
        return false;
      }
    }
  }
  return true;
}

//TODO would be premium if we actually get which configuration is the compatible one
bool PllSplit::compatible(const PllSplit& other) const {
  assert(splitValid());
  assert(other.splitValid());
  return !intersectionSize(other, 0, 1) 
      || !intersectionSize(other, 1, 0) || !intersectionSize(other, 0, 0);
}
//This wonderful monstrosity of an abomination is a quick hack to get which intersection is actually the
//compatible one. 1 for B_1 + A_2; 2 for A_1 + B_2; 3 for B_1 + B_2; 0 for incompatible
int PllSplit::compatiblePREMIUM(const PllSplit& other) const {
  if(!intersectionSize(other, 0, 1)) {return 1;}
  if(!intersectionSize(other, 1, 0)) {return 2;}
  if(!intersectionSize(other, 0, 0)) {return 3;}
  return 0;
}

//TODO ths mask is correct as balls but it should be invoked on the last register only
pll_split_base_t PllSplit::bitmaskForUnusedBits() const {
  pll_split_base_t bit_mask = 0;
  size_t offset = PllSplit::getTipCount() - ((PllSplit::getSplitLen() - 1) * computSplitBaseSize());
  //TODO@ Robin, do one shift instead of a loop
  for(size_t i = 0; i < offset; ++i){
    bit_mask |= (1 << i);
  }
  return bit_mask;
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
bool PllSplit::splitValid() const {
  //This condition sometimes fails on the splits returned from pll lib, needs to be examined!
  //return (_split != nullptr) && !(_split[0] & ~bitmaskForUnusedBits()) && _split[0] & 1u;
  return (_split != nullptr) &&  _split[0] & 1u;
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


//TODO this can be removed if we don't use the other implementation
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
//TODO maybe remove copy and compare references only and access splits directly
bool operator == (const PllSplitList& p1, const PllSplitList& p2) {
  const std::vector<PllSplit>& splits1 = p1.getSplits();
  const std::vector<PllSplit>& splits2 = p2.getSplits();
  if(splits1.size() != splits2.size()) return false;
  for(size_t i = 0; i < splits1.size(); ++i){
    if (!(splits1[i] == splits2[i])) return false;
  }
  return true; 
}
