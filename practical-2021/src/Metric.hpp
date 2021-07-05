#pragma once
#include <string>
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllPointerMap.hpp"
#include "datastructures/IntersectionCache.hpp"
#include "datastructures/TempManager.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"

using PllPosition = size_t;

enum Mode{SIMILARITY, ABSOLUTE, RELATIVE};
static const char *ModeString[] = {"SIMILARITY", "ABSOLUTE", "RELATIVE"};

class Metric {
    public:
    virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, Mode mode, const PllPointerMap& map) const = 0;
    virtual double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const = 0;
    virtual std::string name() const = 0;
    virtual ~Metric() {

    }
};


class GeneralizedMetric : public Metric {
public:
  virtual double evaluate(const PllPosition& s1, const PllPosition& s2, const PllPointerMap& map) const = 0;
  virtual ~GeneralizedMetric() override {

    }

  double distanceOf(const PllSplitList& first, const PllSplitList& second, Mode mode, const PllPointerMap& map) const override {
    std::vector<std::vector<double>> similarities = similaritiesForSplits(first, second, map);
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

    std::vector<std::vector<double>> similaritiesForSplits(const PllSplitList& first, const PllSplitList& second, const PllPointerMap& map) const{
      assert(first.getSplits().size() == first.getSplits().size());
      size_t n = first.getSplits().size();
      std::vector<std::vector<double>>  result = std::vector<std::vector<double>>(n, std::vector<double>(n));
      for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j) {
          //map[first.pos(i)];
          //std::cout << first.pos(i) << first[i]->toString() << "---------------------\n";
          result[i][j] = evaluate(first[i], second[j], map);
          //result[i][j] = evaluate(first[i], second[j]);
        }
      }
      return result;
    }
};

class MSIMetric : public GeneralizedMetric {
  public:
  double evaluate(const PllPosition& s1, const PllPosition& s2, const PllPointerMap& map) const override {
    const PllSplit& sp1 = map[s1];
    if (s1 == s2) return sp1.h();
    const PllSplit& sp2 = map[s2];
    return std::max(phylomath::h(sp1.intersectionSize(sp2, Block_A, Block_A), sp1.intersectionSize(sp2, Block_B, Block_B)),
                    phylomath::h(sp1.intersectionSize(sp2, Block_B, Block_A), sp1.intersectionSize(sp2, Block_A, Block_B)));
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
  double evaluate(const PllPosition& s1, const PllPosition& s2, const PllPointerMap& map) const override {
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    const PllSplit& sp1 = map[s1];
    const PllSplit& sp2 = map[s2];
    assert(sp1.intersectionSize(sp2, Block_A, Block_A) > 0);

    size_t a_1 = sp1.partitionSizeOf(Block_A);
    size_t a_2 = sp2.partitionSizeOf(Block_A);
    size_t b_1 = sp1.partitionSizeOf(Block_B);
    size_t b_2 = sp2.partitionSizeOf(Block_B);

    if (s1 == s2) return phylomath::h(a_2, b_2);

    double phylo_shared;
    size_t intersect_b_a = sp1.intersectionSize(sp2, Block_B, Block_A);
    if(!intersect_b_a) {
      phylo_shared = phylomath::h(b_1, a_2, a_1 + b_1);
    } else {
      size_t intersect_b_b = b_1 - intersect_b_a;
      if(!intersect_b_b){
        phylo_shared = phylomath::h(b_1, b_2, a_1 + b_1);
      } else {
        size_t intersect_a_b = b_2 - intersect_b_b;
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
    double evaluate(const PllPosition& s1, const PllPosition& s2, const PllPointerMap& map) const override {
      const PllSplit& sp1 = map[s1];
      const PllSplit& sp2 = map[s2];
      
      size_t size_a1 = sp1.partitionSizeOf(Block_A);
      size_t size_b1 = sp1.partitionSizeOf(Block_B);
      size_t size_a2 = sp2.partitionSizeOf(Block_A);
      size_t size_b2 = sp2.partitionSizeOf(Block_B);
      return mutualInformation(sp1, Block_A, sp2, Block_A) + mutualInformation(sp1, Block_A, sp2, Block_B) 
           + mutualInformation(sp1, Block_B, sp2, Block_A) + mutualInformation(sp1, Block_B, sp2, Block_B);
    }
    double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
      return plist1.getMaximumEntropy() + plist2.getMaximumEntropy();
    }

    std::string name() const override {
      return "MCI";
    }



    private:
    /*double mutualInformation(const size_t intersection_size, const size_t size_of_partition_block1,
                             const size_t size_of_partition_block2) const {
      assert(size_of_partition_block1 > 0 && size_of_partition_block2 > 0);
      if(intersection_size == 0) {
        return 0.0;
      }
      assert(intersection_size > 0);
      //TODO most of this could be cached
      double pcl = phylomath::clusteringProbability(intersection_size);
      //double p_1 = phylomath::clusteringProbability(size_of_partition_block1);
      //double p_2 = phylomath::clusteringProbability(size_of_partition_block2);
      //double res = pcl * std::log2(pcl / (p_1 * p_2));
      return pcl * phylomath::quickerClusteringProbability(intersection_size, size_of_partition_block1, size_of_partition_block2);
    }*/
    double mutualInformation(const PllSplit& s1,
                             const Partition block_s1, const PllSplit& s2, const Partition block_s2) const {
        //This is a hardcoded statement. The math agrees that x log(x) -> 0 but c++ refuses
        size_t intersection_size = s1.intersectionSize(s2, block_s1, block_s2);
        if(intersection_size == 0){
          return 0.0;
        }
        //TODO most of this could be cached
        double pcl = phylomath::clusteringProbability(intersection_size);
        assert(pcl > 0);
        double p_1 = phylomath::clusteringProbability(s1.partitionSizeOf(block_s1));
        double p_2 = phylomath::clusteringProbability(s2.partitionSizeOf(block_s2));
        assert(p_1 > 0 && p_2 > 0);
        return pcl * std::log2(pcl / (p_1 * p_2));
    }
};


class RFMetric : public Metric {
public:
  virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, Mode mode, const PllPointerMap& map) const override {
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
