#include "Metric.hpp"
#include "DistanceUtil.hpp"
#include "PhylogeneticMathUtils.hpp"
class MSI : public Metrics {
  public: 
  double evaluate(const PllSplit& s1, const PllSplit& s2) const override {
    if (s1 == s2) return phylomath::h(s1); // At this spot we would check for h(A_1, A_2), h(A_1, B_2) and one of these should be 0: Rethink whether dangerous
    return std::max(phylomath::h(s1.intersectionSize(s2, 1, 1), s1.intersectionSize(s2, 0, 0)),
                    phylomath::h(s1.intersectionSize(s2, 0, 1), s1.intersectionSize(s2, 1, 0)));
  }  
  double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
    //TODO can we really expect both lists to be equally long...seems weird
    //TODO assert what the line above critiques (if we decide that splitlists are equally long)
    double result = 0;
    for(unsigned i = 0; i < plist1.getSplitCount(); ++i) {
      result += phylomath::h(plist1[i]);
      result += phylomath::h(plist2[i]);
    }
    return result;
  }
};