#pragma once
#include <string>
#include "datastructures/PllSplits.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"

enum Mode{SIMILARITY, ABSOLUTE, RELATIVE};
static const char *ModeString[] = {"SIMILARITY", "ABSOLUTE", "RELATIVE"};

class Metric {
    public:
    virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, Mode mode) const = 0;
    virtual double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const = 0;
    virtual std::string name() const = 0;
    virtual ~Metric() {

    }
};


class GeneralizedMetric : public Metric {
public:
  virtual double evaluate(const PllSplit* s1, const PllSplit* s2) const = 0;
  virtual ~GeneralizedMetric() override {

    }

  double distanceOf(const PllSplitList& first, const PllSplitList& second, Mode mode) const override {
    std::vector<std::vector<double>> similarities = similaritiesForSplits(first, second);
    assert(similarities.size() == first.getSplits().size());
    double similarity = MaximumMatcher::match(similarities);
    assert(similarity >= 0);
    return (mode == SIMILARITY) ? similarity : distanceFromSimilarity(first, second, similarity, mode);
  }


    double distanceFromSimilarity(const PllSplitList& first,
                                  const PllSplitList& second, double similarity, Mode mode) const{
      double max_value = maximum(first, second);
      assert((max_value - 2 * similarity) > -0.000001);
      double dist = std::max(0.0, max_value - 2 * similarity);
      return (mode == RELATIVE) ? dist : (dist / max_value);
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

class MSIMetric : public GeneralizedMetric {
  public:
  double evaluate(const PllSplit* s1, const PllSplit* s2) const override {
    if (s1 == s2) return s1->h();
    size_t intersect_a_a = s1->intersectionSize(s2);
    size_t intersect_b_a = s2->partitionSizeOf(Block_A) - intersect_a_a;
    size_t intersect_b_b = s1->partitionSizeOf(Block_B) - intersect_b_a;
    size_t intersect_a_b = s2->partitionSizeOf(Block_B) - intersect_b_b;
    return std::max(phylomath::h(intersect_a_a, intersect_b_b),
                    phylomath::h(intersect_b_a, intersect_a_b));
  }
  double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
    return plist1.getMaximumInformationContent() + plist2.getMaximumInformationContent();
  }

  std::string name() const override {
    return "MSI";
  }
};

class SPIMetric : public GeneralizedMetric {
  public:
  double evaluate(const PllSplit* s1, const PllSplit* s2) const override {
    size_t intersect_a_a = s1->intersectionSize(s2);
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(intersect_a_a > 0);

    size_t a_1 = s1->partitionSizeOf(Block_A);
    size_t a_2 = s2->partitionSizeOf(Block_A);
    size_t b_1 = s1->partitionSizeOf(Block_B);
    size_t b_2 = s2->partitionSizeOf(Block_B);

    if (s1 == s2) return phylomath::h(a_2, b_2);

    double phylo_shared;
    size_t intersect_b_a = a_2 - intersect_a_a;
    if(!intersect_b_a) {
      phylo_shared = phylomath::h(b_1, a_2, a_1 + b_1);
    } else {
      size_t intersect_b_b = b_1 - intersect_b_a;
      if(!intersect_b_b){
        phylo_shared = phylomath::h(b_1, b_2, a_1 + b_1);
      } else {
        size_t intersect_a_b = b_2- intersect_b_b;
        if(!intersect_a_b){
          phylo_shared = phylomath::h(a_1, b_2, a_1 + b_1);
        } else {
          //partitions incompatible
          return 0.0;
        }
      }
    }
    return phylomath::h(a_1, b_1) + phylomath::h(a_2, b_2) - phylo_shared;
  }

  double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
    return plist1.getMaximumInformationContent() + plist2.getMaximumInformationContent();
  }

  std::string name() const override {
    return "SPI";
  }

};



class MCIMetric : public GeneralizedMetric {
    public:
    double evaluate(const PllSplit* s1, const PllSplit* s2) const override {
      size_t a_1 = s1->partitionSizeOf(Block_A);
      size_t a_2 = s2->partitionSizeOf(Block_A);
      size_t b_1 = s1->partitionSizeOf(Block_B);
      size_t b_2 = s2->partitionSizeOf(Block_B);
      size_t intersect_a_a = s1->intersectionSize(s2);
      size_t intersect_b_a = a_2 - intersect_a_a;
      size_t intersect_b_b = b_1 - intersect_b_a;
      size_t intersect_a_b = b_2 - intersect_b_b;
        return mutualInformation(a_1, a_2, intersect_a_a) + mutualInformation(b_1, a_2, intersect_b_a)
             + mutualInformation(a_1, b_2, intersect_a_b) + mutualInformation(b_1, b_2, intersect_b_b);
    }
    double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
      return plist1.getMaximumEntropy() + plist2.getMaximumEntropy();
    }

    std::string name() const override {
      return "MCI";
    }



    private:
    double mutualInformation(size_t block_size1, size_t block_size2, size_t intersection_size) const {
        //This is a hardcoded statement. The math agrees that x log(x) -> 0 but c++ refuses
        if(intersection_size == 0){
          return 0.0;
        }
        double pcl = phylomath::clusteringProbability(intersection_size);
        assert(pcl > 0);
        double p_1 = phylomath::clusteringProbability(block_size1);
        double p_2 = phylomath::clusteringProbability(block_size2);
        assert(p_1 > 0 && p_2 > 0);
        return pcl * std::log2(pcl / (p_1 * p_2));
    }
};


class RFMetric : public Metric {
public:
  virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, Mode mode) const override {
    size_t split_count1 = plist1.getSplitCount();
    size_t split_count2 = plist2.getSplitCount();
    if (split_count1 == 0) return static_cast<double> (split_count2);
    if (split_count2 == 0) return static_cast<double> (split_count1);
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
    return (mode == RELATIVE) ? (static_cast<double>(distance) / maximum(plist1, plist2))
                              : static_cast<double>(distance);
  }
  /*OK this is a @Softwipe hack, in order to remove the unused parameter warning I had to remove the name
  but since the signature is still the same its still an override */
  double maximum(const PllSplitList&, const PllSplitList&) const override {
    assert(PllSplit::getTipCount() > 3);
    return static_cast<double> (2 * (PllSplit::getTipCount() - 3));
  }


  virtual std::string name() const override {
    return "RF";
  }

};
