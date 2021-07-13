#include "PllSplits.hpp"
#include "PllTree.hpp"
#include "../PhylogeneticMathUtils.hpp"

size_t PllSplit::tip_count = 0;
size_t PllSplit::split_len = 0;
pll_split_base_t PllSplit::bitmask_for_unused_bits = 0;

PllSplit::PllSplit(pll_split_t s) : _split{s} {
  size_block_A = this->popcount();
  size_block_B = PllSplit::getTipCount() - this->popcount();
  h_value = phylomath::h(size_block_A, size_block_B);
  assert(h_value >= 0);
  entropy_value = phylomath::entropy(size_block_A, size_block_B);
  assert(entropy_value >= 0);
}


bool operator == (const PllSplit& p1, const PllSplit& p2) {
  for(unsigned i = 0; i < PllSplit::split_len; ++i) {
    if(p1._split[i] != p2._split[i]) {
      return false;
    }
  }
  return true;
}
bool operator < (const PllSplit&p1, const PllSplit& p2) {
  for(unsigned i = 0; i < PllSplit::split_len; ++i) {
    if(p1._split[i] != p2._split[i]) {
      return (p1._split[i] < p2._split[i]);
    }
  }
  return false;
}

size_t PllSplit::popcount() const{
  size_t popcount = 0;
  for(size_t i = 0; i < PllSplit::split_len; ++i){
    popcount += __builtin_popcount(_split[i]);
  }
  return popcount;
}

size_t PllSplit::intersectionSize(const PllSplit& other) const {
  pll_split_t other_split = other();
  size_t count = 0;

  for (size_t i = 0; i < PllSplit::split_len - 1; ++i){
    count += __builtin_popcount(_split[i] & other_split[i]);
  }
  count += __builtin_popcount(_split[PllSplit::split_len - 1] & other_split[PllSplit::split_len - 1]
                        & PllSplit::bitmask_for_unused_bits);
  return count;
}

std::string PllSplit::toString() const {
  std::stringstream ss;
  ss << this << ": " << _split << ": ";
  for (size_t i = 0; i < PllSplit::split_len; ++i){
    auto str = std::bitset<32>(_split[i]).to_string();
    std::reverse(str.begin(), str.end());
    ss << str << "|";
  }
  ss << std::endl;
  return ss.str();
}


//TODO find out if preallocation can be done THIS SHOULD NO LONGER BE CALLED
PllSplitList::PllSplitList(const std::vector<size_t> &splits) {
  _split_offsets = splits;
}

bool operator == (const PllSplitList& p1, const PllSplitList& p2) {
  const std::vector<size_t>& splits1 = p1.getSplits();
  const std::vector<size_t>& splits2 = p2.getSplits();
  if(splits1.size() != splits2.size()) return false;
  for(size_t i = 0; i < splits1.size(); ++i){
    if (!(splits1[i] == splits2[i])) return false;
  }
  return true;
}

  void PllSplitList::push(const PllSplit& split, size_t offset) {
    _split_offsets.push_back(offset);
    maximum_entropy += split.entropy();
    maximum_information_content += split.h();
  }

  //Disabled for the restructuring of PLlSplitlists since I wanna remove _splits
/*std::string PllSplitList::toString() const {
  std::stringstream ss;
  ss <<  "-------------------------" << std::endl;
  for(size_t i = 0; i < _splits.size(); ++i){
    ss << _splits[i]->toString();
  }
  ss << "-------------------------" << std::endl;
  return ss.str();
}*/
