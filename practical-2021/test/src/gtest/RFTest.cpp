#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/datastructures/PllPointerMap.hpp"
#include "../../../src/Metric.hpp"
#include "../../../src/datastructures/SimilarityCache.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../../../src/Distances.hpp"
#include "../TestUtil.hpp"

class RFTest : public testing::Test {
  protected:


  std::string current_data_dir = "../test/res/data/";

};






/*TEST_F(RFTest, rf_test) {
  PllSplit::setTipCount(64);
  std::vector<std::vector<size_t>> fst_part1s = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 1, 2, 3, 4, 5},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3}
            };
  PllSplitList fst_splitlist = TestUtil::createSplitList(fst_part1s);

  std::vector<std::vector<size_t>> snd_part1s = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList snd_splitlist = TestUtil::createSplitList(snd_part1s);

  std::vector<std::vector<size_t>> trd_part1s = {
                { 0, 2, 3, 4, 5, 6, 7},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList trd_splitlist = TestUtil::createSplitList(trd_part1s);

  RFMetric rf;
  ASSERT_EQ(rf.distanceOf(fst_splitlist, snd_splitlist, ABSOLUTE), 6);
  ASSERT_EQ(rf.distanceOf(snd_splitlist, fst_splitlist, ABSOLUTE), 6);

  ASSERT_EQ(rf.distanceOf(fst_splitlist, trd_splitlist, ABSOLUTE), 9);
  ASSERT_EQ(rf.distanceOf(trd_splitlist, fst_splitlist, ABSOLUTE), 9);

  ASSERT_EQ(rf.distanceOf(snd_splitlist, trd_splitlist, ABSOLUTE), 3);
  ASSERT_EQ(rf.distanceOf(trd_splitlist, snd_splitlist, ABSOLUTE), 3);

}*/
