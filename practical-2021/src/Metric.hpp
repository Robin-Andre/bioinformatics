#pragma once
#include <string>
#include "datastructures/PllSplits.hpp"
#include "datastructures/UniquePllMap.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"

using PllPosition = size_t;

enum Mode{SIMILARITY, ABSOLUTE, RELATIVE};
static const char *ModeString[] = {"SIMILARITY", "ABSOLUTE", "RELATIVE"};

/**
 * Reprensents a metric for measuring distances amonng phylogenetic trees
 */
class Metric {
    public:
    /**
    * maximum value two trees can admit in this metric
    *
    * @param l1: split list representing the first tree
    * @param l2: split list representing the second tree
    * @return maximum possible value
    */
    virtual double maximum(const PllSplitList& l1, const PllSplitList& l2) const = 0;
    /**
    * @return string identifying this metric
    */
    virtual std::string name() const = 0;
    virtual ~Metric();
};

/**
 * Reprensents a generalized Robinson-Foulds-Distance as presented in https://doi.org/10.1093/bioinformatics/btaa614
 */
class GeneralizedMetric : public Metric {
public:
  /**
  * Evaluate the metric for a pair of splits
  *
  * @param pos1: position fo the first split in the pointer map
  * @param pos2: position fo the second split in the pointer map
  * @param map: pointer map storing the splits
  * @return value for the splits in the metric
  */
  virtual double evaluate(const PllPosition& pos1, const PllPosition& pos2, const UniquePllMap& map) const = 0;
  virtual ~GeneralizedMetric() override;

};

/**
 * Matching Split Information
 * For details refer to https://doi.org/10.1093/bioinformatics/btaa614
 */
class MSIMetric : public GeneralizedMetric {
  public:
  double evaluate(const PllPosition& pos1, const PllPosition& pos2, const UniquePllMap& map) const override {
    const PllSplit& sp1 = map[pos1];
    if (pos1 == pos2) return sp1.h();
    const PllSplit& sp2 = map[pos2];

    size_t intersect_11 = sp1.intersectionSize(sp2);
    size_t intersect_01 = sp2.partitionSizeOf(1) - intersect_11;
    size_t intersect_00 = sp1.partitionSizeOf(0) - intersect_01;
    size_t intersect_10 = sp2.partitionSizeOf(0) - intersect_00;

    assert(intersect_11 <= sp1.partitionSizeOf(1) && intersect_11 <= sp2.partitionSizeOf(1));
    assert(intersect_01 <= sp1.partitionSizeOf(0) && intersect_01 <= sp2.partitionSizeOf(1));
    assert(intersect_00 <= sp1.partitionSizeOf(0) && intersect_00 <= sp2.partitionSizeOf(0));
    assert(intersect_10 <= sp1.partitionSizeOf(1) && intersect_10 <= sp2.partitionSizeOf(0));

    return std::max(phylomath::h(intersect_11, intersect_00),
                    phylomath::h(intersect_01, intersect_10));

  }
  double maximum(const PllSplitList& l1, const PllSplitList& l2) const override {
    return l1.getMaximumInformationContent() + l2.getMaximumInformationContent();
  }

  std::string name() const override {
    return "MSI";
  }

  virtual ~MSIMetric() override;
};


/**
 * Shared Phylogenetic Information
 * For details refer to https://doi.org/10.1093/bioinformatics/btaa614
 */
class SPIMetric : public GeneralizedMetric {
  public:
  double evaluate(const PllPosition& pos1, const PllPosition& pos2, const UniquePllMap& map) const override {
    //because of normalization, the 1-Partitions of pos1 and pos2 always overlap
    const PllSplit& sp1 = map[pos1];
    const PllSplit& sp2 = map[pos2];


    size_t s1_size_1 = sp1.partitionSizeOf(1);
    size_t s2_size_1 = sp2.partitionSizeOf(1);
    size_t s1_size_0 = sp1.partitionSizeOf(0);
    size_t s2_size_0 = sp2.partitionSizeOf(0);

    if (pos1 == pos2) return phylomath::h(s2_size_1, s2_size_0);

    double phylo_shared;
    size_t intersect_11 = sp1.intersectionSize(sp2);
    size_t intersect_01 = s2_size_1 - intersect_11;

    assert(intersect_11 <= s1_size_1 && intersect_11 <= s2_size_1);
    assert(intersect_01 <= s1_size_0 && intersect_01 <= s2_size_1);

    if(!intersect_01) {
      phylo_shared = phylomath::h_shared(s1_size_0, s2_size_1);
    } else {
      size_t intersect_00 = s1_size_0 - intersect_01;
      assert(intersect_00 <= s1_size_0 && intersect_00 <= s2_size_0);
      if(!intersect_00){
        phylo_shared = phylomath::h_shared(s1_size_0, s2_size_0);
      } else {
        size_t intersect_10 = s2_size_0 - intersect_00;
        assert(intersect_10 <= s1_size_1 && intersect_10 <= s2_size_0);
        if(!intersect_10){
          phylo_shared = phylomath::h_shared(s1_size_1, s2_size_0);
        } else {
          //partitions incompatible
          return 0.0;
        }
      }
    }
    return phylomath::h(s1_size_1, s1_size_0) + phylomath::h(s2_size_1, s2_size_0) - phylo_shared;
  }

  double maximum(const PllSplitList& l1, const PllSplitList& l2) const override {
    return l1.getMaximumInformationContent() + l2.getMaximumInformationContent();
  }

  std::string name() const override {
    return "SPI";
  }

  virtual ~SPIMetric() override;

};


/**
 * Mutual Clustering Information
 * For details refer to https://doi.org/10.1093/bioinformatics/btaa614
 */
class MCIMetric : public GeneralizedMetric {
    public:
    double evaluate(const PllPosition& pos1, const PllPosition& pos2, const UniquePllMap& map) const override {
      const PllSplit& sp1 = map[pos1];
      const PllSplit& sp2 = map[pos2];

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
    double maximum(const PllSplitList& l1, const PllSplitList& l2) const override {
      return l1.getMaximumEntropy() + l2.getMaximumEntropy();
    }

    std::string name() const override {
      return "MCI";
    }

    virtual ~MCIMetric() override;



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
      return pcl * phylomath::clusteringProbability(intersection_size, size_of_partition_block1, size_of_partition_block2);
    }

};

/**
  * Common Robinson-Foulds-Distance
  */

class RFMetric : public Metric {
public:
  virtual double distanceOf(const PllSplitList& l1, const PllSplitList& l2, Mode mode) const {
    size_t split_count1 = l1.getSplitCount();
    size_t split_count2 = l2.getSplitCount();
    if (split_count1 == 0) return static_cast<double> (split_count2);
    if (split_count2 == 0) return static_cast<double> (split_count1);
    size_t i = 0;
    size_t j = 0;
    size_t distance = 0;
    while (i < split_count1 && j < split_count2){
      if (l1[i] == l2[j]) {
        ++i;
        ++j;
      } else if (l1[i] < l2[j]) {
        ++distance;
        ++i;
      } else {
        ++distance;
        ++j;
      }
    }
    distance += (split_count1 - i);
    distance += (split_count2 - j);
    return (mode == RELATIVE) ? (static_cast<double>(distance) / maximum(l1, l2))
                              : static_cast<double>(distance);
  }
  /* This override does not name the parameters of the PllSplitLists since they are not required for
     standard RF-Distance */
  double maximum(const PllSplitList&, const PllSplitList&) const override {
    assert(PllSplit::getTipCount() > 3);
    return static_cast<double> (2 * (PllSplit::getTipCount() - 3));
  }


  virtual std::string name() const override {
    return "RF";
  }

};
