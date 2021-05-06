#pragma once
#include "../../src/RFDistance.hpp"
#include "../../src/PllSplits.hpp"

class TestUtil {

public:

  static void rf_data_eq(const RFData& d1, const RFData& d2, float epsilon) {
    EXPECT_NEAR(d1.average_distance, d2.average_distance, epsilon);
    EXPECT_EQ(d1.distances.size(), d2.distances.size());
    EXPECT_EQ(d1.unique_count, d2.unique_count);
    for(unsigned i = 0; i < d1.distances.size(); ++i) {
      EXPECT_EQ(d1.distances[i], d2.distances[i]);
      EXPECT_NEAR(d1.relative_distances[i], d2.relative_distances[i], epsilon);
    }
  }


  static void  split_lists_eq(const PllSplitList& l1, const PllSplitList& l2){
    std::vector<PllSplit> splits1 = l1.getSplits();
    std::vector<PllSplit> splits2 = l2.getSplits();
    EXPECT_EQ(splits1.size(), splits2.size());
    for(size_t i = 0; i < splits1.size(); ++i){
      EXPECT_EQ(splits1[i], splits2[i]);
    }

  }
  static void split_vector_eq(const std::vector<PllSplit>& l1, const std::vector<PllSplit>& l2) {
    EXPECT_EQ(l1.size(), l2.size());
    for(size_t i = 0; i < l1.size(); ++i){
      EXPECT_EQ(l1[i], l2[i]);
    }

  }
};
