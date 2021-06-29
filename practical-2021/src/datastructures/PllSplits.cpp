#include "PllSplits.hpp"
#include "PllTree.hpp"
#include "../enums.hpp"
#include "../PhylogeneticMathUtils.hpp"

size_t PllSplit::tip_count = 0;
size_t PllSplit::split_len = 0;
pll_split_base_t PllSplit::bitmask_for_unused_bits = 0;

PllSplit::PllSplit(pll_split_t s) : _split{s} {
  //assert(splitValid());
  size_block_A = this->popcount();
  size_block_B = PllSplit::getTipCount() - this->popcount();
  h_value = phylomath::h(size_block_A, size_block_B);
  entropy_value = phylomath::entropy(size_block_A, size_block_B);
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
  //assert(splitValid());
  size_t popcount = 0;
  for(size_t i = 0; i < PllSplit::split_len; ++i){
    popcount += __builtin_popcount(_split[i]);
  }
  return popcount;
}
//TODO, we should kill this once we are certain that we don't need it
uint32_t PllSplit::bitExtract(size_t bit_index) const {
  //assert(splitValid());
  assert(bit_index < PllSplit::getTipCount());
  pll_split_base_t split_part = _split[computeMajorIndex(bit_index)];
  return (split_part & (1u << computeMinorIndex(bit_index))) >> computeMinorIndex(bit_index);
}


size_t PllSplit::intersectionSize(const PllSplit& other,
                                  Partition partition_this, Partition partition_other) const {
  //assert(splitValid());
  //assert(other.splitValid());
  pll_split_t other_split = other();
  pll_split_base_t this_mask = (partition_this == Block_A) ? 0 : ~0u; //This is a xor mask so it is flipped for A/B
  pll_split_base_t other_mask = (partition_other == Block_A) ? 0 : ~0u;
  size_t count = 0;

  for (size_t i = 0; i < PllSplit::split_len - 1; ++i){
    count += __builtin_popcount((_split[i] ^ this_mask) & (other_split[i] ^ other_mask));
  }
  count += __builtin_popcount((_split[PllSplit::split_len - 1] ^ this_mask) & (other_split[PllSplit::split_len - 1] ^ other_mask)
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

//TODO find out if preallocation can be done.
PllSplitList::PllSplitList(const PllTree &tree) {
  assert(PllSplit::getTipCount() == tree.getTipCount());
  pll_split_t* tmp_splits = pllmod_utree_split_create(tree.tree()->vroot, tree.getTipCount(), nullptr);
  //_splits.reserve(tree.tree()->tip_count - 3);
  maximum_entropy = 0.0;
  maximum_information_content = 0.0;
  for (size_t i = 0; i < tree.getTipCount() - 3; ++i) {
    _splits.emplace_back(new PllSplit(tmp_splits[i]));
    maximum_entropy += _splits.back()->entropy();
    maximum_information_content += _splits.back()->h();
  }
  free(tmp_splits);

}
//TODO find out if preallocation can be done
PllSplitList::PllSplitList(const std::vector<PllSplit*> &splits) {
  if(splits.size() > 0){
    size_t split_len = PllSplit::getSplitLen();
    pll_split_t split_pointer = static_cast<pll_split_t> (calloc(splits.size()* split_len, sizeof(pll_split_base_t)));
    maximum_entropy = 0.0;
    maximum_information_content = 0.0;
    for (size_t i=0; i<splits.size(); ++i) {
      //TODO rethink if correct
      memcpy(split_pointer + i* split_len, splits[i], split_len * sizeof(pll_split_base_t));
      _splits.emplace_back(new PllSplit(split_pointer + i*split_len));
      maximum_entropy += _splits.back()->entropy();
      maximum_information_content += _splits.back()->h();
    }
  }
}

PllSplitList::~PllSplitList() {
  if (!_splits.empty()) { /*free(_splits[0]);*/ }
}
bool operator == (const PllSplitList& p1, const PllSplitList& p2) {
  const std::vector<PllSplit*>& splits1 = p1.getSplits();
  const std::vector<PllSplit*>& splits2 = p2.getSplits();
  if(splits1.size() != splits2.size()) return false;
  for(size_t i = 0; i < splits1.size(); ++i){
    if (!(splits1[i] == splits2[i])) return false;
  }
  return true;
}

  void PllSplitList::push(PllSplit* split) {
    _splits.push_back(split);
    maximum_entropy += split->entropy();
    //std::cout << "Push command: " << split << "\n";
    //std::cout << "Dereferenced: " << (*split).toString();
    maximum_information_content += split->h();
    //std::cout << "after: " << split << "\n";
  }
std::string PllSplitList::toString() const {
  std::stringstream ss;
  ss <<  "-------------------------" << std::endl;
  for(size_t i = 0; i < _splits.size(); ++i){
    ss << _splits[i]->toString();
  }
  ss << "-------------------------" << std::endl;
  return ss.str();
}
