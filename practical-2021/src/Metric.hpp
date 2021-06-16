#pragma once
#include <string>
#include "datastructures/PllSplits.hpp"
#include "PhylogeneticMathUtils.hpp"
#include "MaximumMatcher.hpp"
#include <mpfr.h>

enum Mode{SIMILARITY, ABSOLUTE, RELATIVE};
static const char *ModeString[] = {"SIMILARITY", "ABSOLUTE", "RELATIVE"};

class Metric {
    public:
    virtual void distanceOf(mpfr_t result, const PllSplitList& plist1, const PllSplitList& plist2, Mode mode) const = 0;
    virtual void maximum(mpfr_t result, const PllSplitList& plist1, const PllSplitList& plist2) const = 0;
    virtual std::string name() const = 0;
    virtual ~Metric() {

    }
};


class GeneralizedMetric : public Metric{
public:
  virtual void evaluate(mpfr_t result, const PllSplit& s1, const PllSplit& s2) const = 0;
  virtual ~GeneralizedMetric() override {

    }

  void distanceOf(mpfr_t result, const PllSplitList& first, const PllSplitList& second, Mode mode) const override {
    assert(first.getSplits().size() == first.getSplits().size());
    size_t n = first.getSplits().size();
    std::vector<std::vector<double>> similarities = similaritiesForSplits(first, second);
    assert(similarities.size() == first.getSplits().size());
    std::vector<size_t> matching = MaximumMatcher::match_vector(similarities);
    mpfr_t temp;
    mpfr_init_set_ui(temp, 0, RND);
    mpfr_set_ui(result, 0, RND);
    for(size_t i = 0; i < matching.size(); ++i) {
      evaluate(temp, first[i], second[matching[i]]);
      mpfr_add(result, result, temp, RND);
    }
    mpfr_clear(temp);
    assert(mpfr_sgn(result) >= 0);
    if (mode != SIMILARITY) {
      distanceFromSimilarity(result, first, second, mode);
    }
  }



    void distanceFromSimilarity(mpfr_t result, const PllSplitList& first,
                                  const PllSplitList& second, Mode mode) const{
      mpfr_t max_val;
      mpfr_init_set_ui(max_val, 0, RND);
      maximum(max_val, first, second);
      //assert((max_value - 2 * similarity) > -0.000001);
      mpfr_mul_ui(result, result, 2.0, RND);
      mpfr_sub(result, max_val, result, RND);
      if (mpfr_sgn(result) < 0){
        mpfr_set_ui(result, 0, RND);
      }
      //TODO: Set at least 0 ?!
      //assert(mpfr_sgn(result) > 0);
      if (mode == RELATIVE) {
        mpfr_div(result, result, max_val, RND);
      }
      mpfr_clear(max_val);
    }

    std::vector<std::vector<double>> similaritiesForSplits(const PllSplitList& first, const PllSplitList& second) const {
      size_t n = first.getSplits().size();
      std::vector<std::vector<double>> similarities = std::vector<std::vector<double>>(n, std::vector<double>(n));
      mpfr_t temp;
      mpfr_init_set_ui(temp, 0, RND);
      for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
          evaluate(temp, first[i], second[j]);
          similarities[i][j] = mpfr_get_d(temp, RND);
        }
      }
      mpfr_clear(temp);
      return similarities;
    }




};

class MSIMetric : public GeneralizedMetric {
  public:
  void evaluate(mpfr_t result, const PllSplit& s1, const PllSplit& s2) const override {
    mpfr_set_ui(result, 0, RND);
    if (s1 == s2) {
      phylomath::h(result, s1);
      return;
    }
    mpfr_t phylo_1, phylo_2;
    mpfr_init_set_ui(phylo_1, 0, RND);
    mpfr_init_set_ui(phylo_2, 0, RND);
    phylomath::h(phylo_1, s1.intersectionSize(s2, Block_A, Block_A), s1.intersectionSize(s2, Block_B, Block_B));
    phylomath::h(phylo_2, s1.intersectionSize(s2, Block_B, Block_A), s1.intersectionSize(s2, Block_A, Block_B));
    mpfr_max(result, phylo_1, phylo_2, RND);
    mpfr_clear(phylo_1);
    mpfr_clear(phylo_2);
  }
  void maximum(mpfr_t result, const PllSplitList& plist1, const PllSplitList& plist2) const override {
    mpfr_set_ui(result, 0, RND);
    mpfr_t temp;
    mpfr_init_set_ui(temp, 0, RND);
    for(unsigned i = 0; i < plist1.getSplitCount(); ++i) {
      phylomath::h(temp, plist1[i]);
      mpfr_add(result, result, temp, RND);
    }
    for(unsigned i = 0; i < plist2.getSplitCount(); ++i) {
      phylomath::h(temp, plist2[i]);
      mpfr_add(result, result, temp, RND);
    }
    mpfr_clear(temp);
  }

  std::string name() const override {
    return "MSI";
  }
};

class SPIMetric : public GeneralizedMetric {
  public:
  void evaluate(mpfr_t result, const PllSplit& s1, const PllSplit& s2) const override {
    //because of normalization, the 1-Partitions of s1 and s2 always overlap
    assert(s1.intersectionSize(s2, Block_A, Block_A) > 0);

    size_t a_1 = s1.partitionSizeOf(Block_A);
    size_t a_2 = s2.partitionSizeOf(Block_A);
    size_t b_1 = s1.partitionSizeOf(Block_B);
    size_t b_2 = s2.partitionSizeOf(Block_B);

    if (s1 == s2) {
      phylomath::h(result, a_2, b_2);
      return;
    }
    mpfr_t phylo_1, phylo_2;
    mpfr_init_set_ui(phylo_1, 0, RND);
    mpfr_init_set_ui(phylo_2, 0, RND);
    if(!s1.intersectionSize(s2, Block_B, Block_A)) {
      phylomath::h(result, b_1, a_2, a_1 + b_1);
    } else if (!s1.intersectionSize(s2, Block_A, Block_B)) {
      phylomath::h(result, a_1, b_2, a_1 + b_1);
    } else if(!s1.intersectionSize(s2, Block_B, Block_B)) {
      phylomath::h(result, b_1, b_2, a_1 + b_1);
    } else {
      mpfr_set_ui(result, 0, RND);
      return;
    }
    phylomath::h(phylo_1, a_1, b_1);
    phylomath::h(phylo_2, a_2, b_2);
    mpfr_add(phylo_1, phylo_1, phylo_2, RND);
    mpfr_sub(result, phylo_1, result, RND);
    mpfr_clear(phylo_1);
    mpfr_clear(phylo_2);

  }

  void maximum(mpfr_t result, const PllSplitList& plist1, const PllSplitList& plist2) const override {
    mpfr_set_ui(result, 0, RND);
    mpfr_t temp;
    mpfr_init_set_ui(temp, 0, RND);
    for(unsigned i = 0; i < plist1.getSplitCount(); ++i) {
      phylomath::h(temp, plist1[i]);
      mpfr_add(result, result, temp, RND);
    }
    for(unsigned i = 0; i < plist2.getSplitCount(); ++i) {
      phylomath::h(temp, plist2[i]);
      mpfr_add(result, result, temp, RND);
    }
    mpfr_clear(temp);
  }

  std::string name() const override {
    return "SPI";
  }

};



class MCIMetric : public GeneralizedMetric {
    public:
    void evaluate(mpfr_t result, const PllSplit& s1, const PllSplit& s2) const override {
      mpfr_t temp;
      mpfr_init_set_ui(temp, 0, RND);
      mpfr_set_ui(result, 0, RND);
      mutualInformation(temp, s1, Block_A, s2, Block_A);
      mpfr_add(result, result, temp, RND);
      mutualInformation(temp, s1, Block_B, s2, Block_A);
      mpfr_add(result, result, temp, RND);
      mutualInformation(temp, s1, Block_A, s2, Block_B);
      mpfr_add(result, result, temp, RND);
      mutualInformation(temp, s1, Block_B, s2, Block_B);
      mpfr_add(result, result, temp, RND);
      mpfr_clear(temp);
    }
    void maximum(mpfr_t result, const PllSplitList& plist1, const PllSplitList& plist2) const override {
      mpfr_t temp;
      mpfr_init_set_ui(temp, 0, RND);
      mpfr_set_ui(result, 0, RND);
      //TODO it feels weird to A) Recalculate this every time. B) not being able to call it on a single splitlist
      for(size_t i = 0; i < plist1.getSplitCount(); ++i){
        phylomath::entropy(temp, plist1[i]);
        mpfr_add(result, result, temp, RND);
      }
      for(size_t i = 0; i < plist2.getSplitCount(); ++i){
        phylomath::entropy(temp, plist2[i]);
        mpfr_add(result, result, temp, RND);
      }
      mpfr_clear(temp);
    }

    std::string name() const override {
      return "MCI";
    }



    private:
    void mutualInformation(mpfr_t result, const PllSplit&s1,
                             const Partition block_s1, const PllSplit& s2, const Partition block_s2) const {
        //This is a hardcoded statement. The math agrees that x log(x) -> 0 but c++ refuses
        if(s1.intersectionSize(s2, block_s1, block_s2) == 0){
          mpfr_set_ui(result, 0, RND);
          return;
        }
        mpfr_t pcl, p_1, p_2;
        mpfr_init_set_ui(pcl, 0, RND);
        mpfr_init_set_ui(p_1, 0, RND);
        mpfr_init_set_ui(p_2, 0, RND);
        phylomath::clusteringProbability(pcl, s1, block_s1, s2, block_s2);
        assert(mpfr_sgn(pcl) > 0);
        phylomath::clusteringProbability(p_1, s1, block_s1);
        phylomath::clusteringProbability(p_2, s2, block_s2);
        assert(mpfr_sgn(p_1) > 0 && mpfr_sgn(p_2) > 0);
        mpfr_mul(p_1, p_1, p_2, RND);
        mpfr_div(p_1, pcl, p_1, RND);
        mpfr_log2(p_1, p_1, RND);
        mpfr_mul(result, pcl, p_1, RND);
        mpfr_clear(p_1);
        mpfr_clear(p_2);
        mpfr_clear(pcl);
    }
};


class RFMetric : public Metric {
public:
  virtual void distanceOf(mpfr_t result, const PllSplitList& plist1, const PllSplitList& plist2, Mode mode) const override {
    size_t split_count1 = plist1.getSplitCount();
    size_t split_count2 = plist2.getSplitCount();
    if (split_count1 == 0) {
      mpfr_set_ui(result, split_count2, RND);
      return;
    }
    if (split_count2 == 0) {
      mpfr_set_ui(result, split_count1, RND);
      return;
    }
    size_t i = 0;
    size_t j = 0;
    mpfr_set_ui(result, 0, RND);
    while (i < split_count1 && j < split_count2){
      if (plist1[i] == plist2[j]) {
        ++i;
        ++j;
      } else if (plist1[i] < plist2[j]) {
        mpfr_add_ui(result, result, 1, RND);
        ++i;
      } else {
        mpfr_add_ui(result, result, 1, RND);
        ++j;
      }
    }
    mpfr_add_ui(result, result, (split_count1 - i), RND);
    mpfr_add_ui(result, result, (split_count2 - j), RND);
    if(mode == RELATIVE) {
      mpfr_t max_val;
      mpfr_init_set_ui(max_val, 0, RND);
      maximum(max_val, plist1, plist2);
      mpfr_div(result, result, max_val, RND);
      mpfr_clear(max_val);
    }

  }
  /*OK this is a @Softwipe hack, in order to remove the unused parameter warning I had to remove the name
  but since the signature is still the same its still an override */
  void maximum(mpfr_t result, const PllSplitList&, const PllSplitList&) const override {
    assert(PllSplit::getTipCount() > 3);
    mpfr_set_ui(result, PllSplit::getTipCount(), RND);
    mpfr_sub_ui(result, result, 3, RND);
    mpfr_mul_ui(result, result, 2, RND);
  }


  virtual std::string name() const override {
    return "RF";
  }

};
