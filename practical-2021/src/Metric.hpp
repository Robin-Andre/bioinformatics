#pragma once
#include <string>
#include "datastructures/PllSplits.hpp"
#include "datastructures/PllPointerMap.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"

using PllPosition = size_t;

enum Mode{SIMILARITY, ABSOLUTE, RELATIVE};
static const char *ModeString[] = {"SIMILARITY", "ABSOLUTE", "RELATIVE"};

class Metric {
    public:
    virtual double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const = 0;
    virtual std::string name() const = 0;
    virtual ~Metric() {}
};


class GeneralizedMetric : public Metric {
public:
  virtual double evaluate(const PllPosition& s1, const PllPosition& s2, const PllPointerMap& map) const = 0;
  virtual ~GeneralizedMetric() override {}

};

class MSIMetric : public GeneralizedMetric {
  public:
  double evaluate(const PllPosition& s1, const PllPosition& s2, const PllPointerMap& map) const override {
    const PllSplit& sp1 = map[s1];
    if (s1 == s2) return sp1.h();
    const PllSplit& sp2 = map[s2];
    size_t intersect_11 = sp1.intersectionSize(sp2);
    assert(intersect_11 <= sp1.partitionSizeOf(1) && intersect_11 <= sp2.partitionSizeOf(1));
    size_t intersect_01 = sp2.partitionSizeOf(1) - intersect_11;
    assert(intersect_01 <= sp1.partitionSizeOf(0) && intersect_01 <= sp2.partitionSizeOf(1));
    size_t intersect_00 = sp1.partitionSizeOf(0) - intersect_01;
    assert(intersect_00 <= sp1.partitionSizeOf(0) && intersect_00 <= sp2.partitionSizeOf(0));
    size_t intersect_10 = sp2.partitionSizeOf(0) - intersect_00;
    assert(intersect_10 <= sp1.partitionSizeOf(1) && intersect_10 <= sp2.partitionSizeOf(0));
    return std::max(phylomath::h(intersect_11, intersect_00),
                    phylomath::h(intersect_01, intersect_10));

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


    size_t s1_size_1 = sp1.partitionSizeOf(1);
    size_t s2_size_1 = sp2.partitionSizeOf(1);
    size_t s1_size_0 = sp1.partitionSizeOf(0);
    size_t s2_size_0 = sp2.partitionSizeOf(0);

    if (s1 == s2) return phylomath::h(s2_size_1, s2_size_0);

    double phylo_shared;
    size_t intersect_11 = sp1.intersectionSize(sp2);
    assert(intersect_11 <= s1_size_1 && intersect_11 <= s2_size_1);
    size_t intersect_01 = s2_size_1 - intersect_11;
    assert(intersect_01 <= s1_size_0 && intersect_01 <= s2_size_1);
    if(!intersect_01) {
      phylo_shared = phylomath::h(s1_size_0, s2_size_1, s1_size_1 + s1_size_0);
    } else {
      size_t intersect_00 = s1_size_0 - intersect_01;
      assert(intersect_00 <= s1_size_0 && intersect_00 <= s2_size_0);
      if(!intersect_00){
        phylo_shared = phylomath::h(s1_size_0, s2_size_0, s1_size_1 + s1_size_0);
      } else {
        size_t intersect_10 = s2_size_0 - intersect_00;
        assert(intersect_10 <= s1_size_1 && intersect_10 <= s2_size_0);
        if(!intersect_10){
          phylo_shared = phylomath::h(s1_size_1, s2_size_0, s1_size_1 + s1_size_0);
        } else {
          //partitions incompatible
          return 0.0;
        }
      }
    }
    return phylomath::h(s1_size_1, s1_size_0) + phylomath::h(s2_size_1, s2_size_0) - phylo_shared;
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

      size_t s1_size_1 = sp1.partitionSizeOf(1);
      size_t s1_size_0 = sp1.partitionSizeOf(0);
      size_t s2_size_1 = sp2.partitionSizeOf(1);
      size_t s2_size_0 = sp2.partitionSizeOf(0);
      size_t intersect_11 = sp1.intersectionSize(sp2);
      assert(intersect_11 <= s1_size_1 && intersect_11 <= s2_size_1);
      size_t intersect_01 = s2_size_1 - intersect_11;
      assert(intersect_01 <= s1_size_0 && intersect_01 <= s2_size_1);
      size_t intersect_00 = s1_size_0 - intersect_01;
      assert(intersect_00 <= s1_size_0 && intersect_00 <= s2_size_0);
      size_t intersect_10 = s2_size_0 - intersect_00;
      assert(intersect_10 <= s1_size_1 && intersect_10 <= s2_size_0);
      return mutualInformation(intersect_11, s1_size_1, s2_size_1) +
             mutualInformation(intersect_10, s1_size_1, s2_size_0) +
             mutualInformation(intersect_01, s1_size_0, s2_size_1) +
             mutualInformation(intersect_00, s1_size_0, s2_size_0);

    }
    double maximum(const PllSplitList& plist1, const PllSplitList& plist2) const override {
      return plist1.getMaximumEntropy() + plist2.getMaximumEntropy();
    }

    std::string name() const override {
      return "MCI";
    }



    private:
    double mutualInformation(const size_t intersection_size, const size_t size_of_partition_block1,
                             const size_t size_of_partition_block2) const {
      assert(size_of_partition_block1 > 0 && size_of_partition_block2 > 0);
      if(intersection_size == 0) {
        return 0.0;
      }
      assert(intersection_size > 0);
      //TODO most of this could be cached
      double pcl = phylomath::clusteringProbability(intersection_size);
      return pcl * phylomath::quickerClusteringProbability(intersection_size, size_of_partition_block1, size_of_partition_block2);
    }

};


class RFMetric : public Metric {
public:
  virtual double distanceOf(const PllSplitList& plist1, const PllSplitList& plist2, Mode mode) const {
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
