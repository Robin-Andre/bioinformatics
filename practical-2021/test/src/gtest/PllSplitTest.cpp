#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../TestUtil.hpp"


class PllSplitTest : public testing::Test {

protected:

  void split_vector_eq(const std::vector<PllSplit*>& l1, const std::vector<PllSplit*>& l2) {
    EXPECT_EQ(l1.size(), l2.size());
    for(size_t i = 0; i < l1.size(); ++i){
      EXPECT_EQ(*l1[i], *l2[i]);
    }
  }
};



TEST_F(PllSplitTest, test_popcount) {
  PllSplit::setTipCount(64);
  // tips in the partition 1 (0 needs to be in)
  std::vector<size_t> part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_EQ(part1.size(), split.partitionSizeOf(1));
  free(split());
}


TEST_F(PllSplitTest, test_operators) {
  PllSplit::setTipCount(64);
  // tips in the partition 1
  // contains 0 impicitly
  std::vector<size_t> fst_part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit fst_split = TestUtil::createSplit(fst_part1);
  std::vector<size_t> snd_part1 = {0, 1, 4, 8, 19, 45, 63};
  PllSplit snd_split = TestUtil::createSplit(snd_part1);
  ASSERT_TRUE(snd_split < fst_split);

  std::vector<size_t> trd_part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit trd_split = TestUtil::createSplit(trd_part1);
  ASSERT_TRUE(trd_split == fst_split);
  ASSERT_TRUE(snd_split < fst_split);
  free(fst_split());
  free(snd_split());
  free(trd_split());
}


TEST_F(PllSplitTest, test_intersectcount) {
    PllSplit::setTipCount(10);
    std::vector<size_t> part1_a = {0, 2, 4, 9};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 1, 7, 9};
    PllSplit split_b = TestUtil::createSplit(part1_b);

    EXPECT_EQ(split_a.intersectionSize(split_a), 4);
    EXPECT_EQ(split_a.intersectionSize(split_b), 2);

    free(split_a());
    free(split_b());
}


//These test cases were added because of the deprecated method of bitmask
TEST_F(PllSplitTest, popcount_2registers_nonfull) {
  PllSplit::setTipCount(36);
  PllSplit test_split = TestUtil::createSplit({0, 31, 33, 35});
  EXPECT_EQ(test_split.partitionSizeOf(1), 4);
  free(test_split());
}
TEST_F(PllSplitTest, intersectionsize_2registers_nonfull) {
  PllSplit::setTipCount(36);
  PllSplit test_split = TestUtil::createSplit({0, 31, 33, 35});
  EXPECT_EQ(test_split.intersectionSize(test_split), 4);
  free(test_split());
}
//Remember that the splits are sorted by LSB, but our print method is reverting the order
TEST_F(PllSplitTest, lessthanoperator_2registers_nonfull) {
  PllSplit::setTipCount(36);
  PllSplit test_split1 = TestUtil::createSplit({0, 31, 33, 35});
  // 1000 0000 0000 0000 0000 0000 0000 0001 | 1010
  PllSplit test_split2 = TestUtil::createSplit({0, 30, 33, 35});
  // 0100 0000 0000 0000 0000 0000 0000 0001 | 1010
  PllSplit test_split3 = TestUtil::createSplit({0, 31, 34, 35});
  // 1000 0000 0000 0000 0000 0000 0000 0001 | 1100
  EXPECT_TRUE(test_split2 < test_split1);
  EXPECT_TRUE(test_split1 < test_split3);
  EXPECT_TRUE(test_split2 < test_split3);
  free(test_split1());
  free(test_split2());
  free(test_split3());
}
TEST_F(PllSplitTest, equaloperator_2registers_nonfull) {
  PllSplit::setTipCount(36);
  PllSplit test_split = TestUtil::createSplit({0, 31, 33, 35});
  EXPECT_TRUE(test_split == test_split);
  free(test_split());
}
