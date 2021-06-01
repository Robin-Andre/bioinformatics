#include "Metric.hpp"
#include "../enums.hpp"
#include "../PhylogeneticMathUtils.hpp"
class MCI : public Metrics {
    public:
    double evaluate(const PllSplit& s1, const PllSplit& s2) const override {
        return helper(s1, Block_A, s2, Block_A) + helper(s1, Block_B, s2, Block_A) 
             + helper(s1, Block_A, s2, Block_B) + helper(s1, Block_B, s2, Block_B);
    }
    double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
      double result;  
      //TODO it feels weird to A) Recalculate this every time. B) not being able to call it on a single splitlist
      for(size_t i = 0; i < plist1.getSplitCount(); ++i){
        result += phylomath::entropy(plist1[i]);
        result += phylomath::entropy(plist2[i]);
      }
      return result;
    }
    private: 
    double helper(const PllSplit&s1, const Partition block_s1, const PllSplit& s2, const Partition block_s2) const {
        double pcl = phylomath::clusteringProbability(s1, block_s1, s2, block_s2);
        //This is a hardcoded statement. The math agrees that x log(x) -> 0 but c++ refuses
        if(pcl == 0) {
            return 0;
        }
        double p_1 = phylomath::clusteringProbability(s1, block_s1);
        double p_2 = phylomath::clusteringProbability(s2, block_s2);
        return pcl * std::log2(pcl / (p_1 * p_2));
    }
};