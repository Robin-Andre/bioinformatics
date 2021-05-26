#include "Metric.hpp"
#include "PhylogeneticMathUtils.hpp"
class SPI : public Metrics {
  public:
  double evaluate(const PllSplit& s1, const PllSplit& s2) const override {
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(s1.intersectionSize(s2, 1, 1) > 0);
    if (!s1.compatible(s2)) return 0;
    size_t a_1 = s1.partitionSize(1);
    size_t a_2 = s2.partitionSize(1);
    size_t b_1 = s1.partitionSize(0);
    size_t b_2 = s2.partitionSize(0);
    if (s1 == s2) return phylomath::h(a_2, b_2);
    return phylomath::h(a_1, b_1) + phylomath::h(a_2, b_2) - phylomath::h(a_1, b_1, a_2, b_2);

  
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