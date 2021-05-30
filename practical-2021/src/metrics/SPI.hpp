#include "Metric.hpp"
#include "../PhylogeneticMathUtils.hpp"
class SPI : public Metrics {
  public:
  //This is the old implementation but until I see definite proof that we use the new one it will stay as 
  // copy paste tool
  /*double evaluate(const PllSplit& s1, const PllSplit& s2) const override {
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(s1.intersectionSize(s2, 1, 1) > 0);
    if (!s1.compatible(s2)) return 0;
    size_t a_1 = s1.partitionSizeOf(Block_A);
    size_t a_2 = s2.partitionSizeOf(Block_A);
    size_t b_1 = s1.partitionSizeOf(Block_B);
    size_t b_2 = s2.partitionSizeOf(Block_B);
    if (s1 == s2) return phylomath::h(a_2, b_2);
    return phylomath::h(a_1, b_1) + phylomath::h(a_2, b_2) - phylomath::h(a_1, b_1, a_2, b_2);

  
  }*/
  double evaluate(const PllSplit& s1, const PllSplit& s2) const override {
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(s1.intersectionSize(s2, 1, 1) > 0);
    int compatible_numeric = s1.compatiblePREMIUM(s2);
    if(compatible_numeric == 0) {
      return 0;
    }
    size_t a_1 = s1.partitionSizeOf(Block_A);
    size_t a_2 = s2.partitionSizeOf(Block_A);
    size_t b_1 = s1.partitionSizeOf(Block_B);
    size_t b_2 = s2.partitionSizeOf(Block_B);
    if (s1 == s2) return phylomath::h(a_2, b_2);
    double phylo_shared;
    switch(compatible_numeric) {
      case 1: 
        phylo_shared = phylomath::h(b_1, a_2, a_1 + b_1);
        break;
      case 2:
        phylo_shared = phylomath::h(a_1, b_2, a_1 + b_1);
        break;
      case 3:
        phylo_shared = phylomath::h(b_1, b_2, a_1 + b_1);
        break;
    }
  
    return phylomath::h(a_1, b_1) + phylomath::h(a_2, b_2) - phylo_shared;

  
  }
  double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
    double result = 0;
    for(unsigned i = 0; i < plist1.getSplitCount(); ++i) {
      result += phylomath::h(plist1[i]);
      result += phylomath::h(plist2[i]);
    }
    return result;
  }
};