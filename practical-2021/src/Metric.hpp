#pragma once
#include <string>
#include "datastructures/PllSplits.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"

class Metric {
    public:
    virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, bool normalize) const = 0;
    virtual std::string name() const = 0;
};


class GeneralizedMetric : public Metric{
public:
  virtual double evaluate(const PllSplit& s1, const PllSplit& s2) const = 0;
  virtual double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const = 0;

  double distanceOf(const PllSplitList& first, const PllSplitList& second, bool normalize) const {
    std::vector<std::vector<double>> similarities = similaritiesForSplits(first, second);
    double similarity = MaximumMatcher::match(similarities);
    /*for (size_t k = 0; k < similarities.size(); ++k){
      for (size_t l = 0; l < similarities[k].size(); ++l){
        std::cout << similarities[k][l] << "; ";
      }
      std::cout << "|" << phylomath::h(tree_splits[i][k]);
      std::cout << std::endl;
    }
    std::cout << "SIM: " << similarity << std::endl;*/
    return normalize ? distanceFromSimilarity(first, second, similarity) : similarity;
  }

    double maximumValue(const PllSplitList& first, const PllSplitList& second) const{
      return (maximum(first, second) / 2);
    }

    double distanceFromSimilarity(const PllSplitList& first, const PllSplitList& second, double similarity) const{
      return maximumValue(first, second) - similarity;
    }

    std::vector<std::vector<double>> similaritiesForSplits(const PllSplitList& first, const PllSplitList& second) const{
      assert(first.getSplits().size() == first.getSplits().size());
      size_t n = first.getSplits().size();
      std::vector<std::vector<double>>  result = std::vector<std::vector<double>>(n, std::vector<double>(n));
      for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
          result[i][j] = evaluate(first[i], second[j]);
        }
      }
      return result;
    }
};

class MSI : public GeneralizedMetric {
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

  std::string name() const override {
    return "MSI";
  }
};

class SPI : public GeneralizedMetric {
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

  std::string name() const override {
    return "SPI";
  }
};



class MCI : public GeneralizedMetric {
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

    std::string name() const override {
      return "MCI";
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


class RF : public Metric {
public:
  virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, bool normalize) const override {
    size_t split_count1 = plist1.getSplitCount();
    size_t split_count2 = plist2.getSplitCount();
    if (split_count1 == 0) return split_count2;
    if (split_count2 == 0) return split_count1;
    size_t i = 0;
    size_t j = 0;
    size_t distance = 0;
    while (i < split_count1 && j < split_count2){
      if (plist1[i] == plist2[j]) {
        ++i;
        ++j;
      } else if (plist1[i] < plist2[j]) {
        ++distance;
        ++i;
      } else {
        ++distance;
        ++j;
      }
    }
    distance += (split_count1 - i);
    distance += (split_count2 - j);
    return distance;
  }


  virtual std::string name() const override {
    return "RF";
  }
};
