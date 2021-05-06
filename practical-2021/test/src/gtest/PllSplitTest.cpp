#include "gtest/gtest.h"
#include "../../../src/PllSplits.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../TestUtil.hpp"


class PllSplitTest : public testing::Test {
protected:
  PllSplit createSplit(std::vector<size_t> part1) {
    if (part1[0] != 0) throw "In every split, 0 must be in Partition 1, hence it must hold that part1[0]==0";
    auto split_bits = (pll_split_t)calloc(PllSplit::getSplitLen(), sizeof(pll_split_base_t));
    setBits(split_bits, part1);
    return PllSplit(split_bits);
  }

  void setBits(pll_split_t split_bits, std::vector<size_t> part1) {
    size_t major_idx;
    size_t minor_idx;
    for(size_t tip : part1){
      major_idx = tip / (sizeof(pll_split_base_t) * 8);
      minor_idx = tip % (sizeof(pll_split_base_t) * 8);
      split_bits[major_idx] |= (1 << minor_idx);
    }
  }

  PllSplitList createSplitList(std::vector<std::vector<size_t>> part1s){
    std::vector<PllSplit> splits;
    pll_split_t split_pointer = (pll_split_t) calloc(part1s.size()* PllSplit::getSplitLen(), sizeof(pll_split_base_t));
    for (size_t i=0; i<part1s.size(); ++i) {
      setBits(split_pointer + i*PllSplit::getSplitLen(), part1s[i]);
      splits.emplace_back(PllSplit(split_pointer + i*PllSplit::getSplitLen()));
    }
    std::sort(splits.begin(), splits.end());
    return PllSplitList(splits);
  }
};


TEST_F(PllSplitTest, test_bitExtract) {
    PllSplit::setSplitLen(2);
    // tips in the partition 1 (0 needs to be in)
    std::vector<size_t> part1 = {0, 2, 4, 8, 19, 45, 63};
    PllSplit split = createSplit(part1);
    std::vector<bool> split_representation(64);
    for(size_t tip : part1) {
      split_representation[tip] = true;
    }
    for(size_t i = 0; i < 64; i++){
      EXPECT_EQ(split_representation[i], split.bitExtract(i));
    }
    free(split());
}


TEST_F(PllSplitTest, test_popcount) {
  PllSplit::setSplitLen(2);
  // tips in the partition 1 (0 needs to be in)
  std::vector<size_t> part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit split = createSplit(part1);
  EXPECT_EQ(part1.size(), split.popcount());
  free(split());
}


TEST_F(PllSplitTest, test_operators) {
  PllSplit::setSplitLen(2);
  // tips in the partition 1
  // contains 0 impicitly
  std::vector<size_t> fst_part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit fst_split = createSplit(fst_part1);
  std::vector<size_t> snd_part1 = {0, 1, 4, 8, 19, 45, 63};
  PllSplit snd_split = createSplit(snd_part1);
  ASSERT_TRUE(snd_split < fst_split);

  std::vector<size_t> trd_part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit trd_split = createSplit(trd_part1);
  ASSERT_TRUE(trd_split == fst_split);
  ASSERT_TRUE(snd_split < fst_split);
  free(fst_split());
  free(snd_split());
  free(trd_split());
}

TEST_F(PllSplitTest, test_list_constructor) {
  PllSplit::setSplitLen(2);
  std::vector<PllSplit> splits;
  std::vector<size_t> fst_part1 = {0, 2, 4, 8, 19, 45, 63};
  splits.emplace_back(createSplit(fst_part1));
  std::vector<size_t> snd_part1 = {0, 1, 4, 8, 19, 45, 63};
  splits.emplace_back(createSplit(snd_part1));
  std::vector<size_t> trd_part1 = {0, 19, 29, 39};
  splits.emplace_back(createSplit(trd_part1));
  PllSplitList split_list = PllSplitList(splits);
  TestUtil::split_vector_eq(splits, split_list.getSplits());
  for(PllSplit split : splits){
    free(split());
  }

}

TEST_F(PllSplitTest, test_tree_constructor) {
  PllTree test_tree = PllTree("((a1, a2), (b1,b2), (c, (d1, d2)));");
  PllSplit::setSplitLen(PllSplit::computeSplitLen(test_tree.getTipCount()));
  PllSplitList split_list_from_tree = PllSplitList(test_tree);
  std::vector<std::vector<size_t>> expected = {
                { 0, 1},
                { 0, 1, 2, 3},
                { 0, 1, 4, 5, 6},
                { 0, 1, 2, 3, 4},
  };
  PllSplitList expected_splitlist = createSplitList(expected);
  TestUtil::split_lists_eq(expected_splitlist, split_list_from_tree);

}


TEST_F(PllSplitTest, test_difference) {
  PllSplit::setSplitLen(2);
  // tips in the partition 1
  // contains 0 impicitly
  std::vector<std::vector<size_t>> fst_part1s = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 1, 2, 3, 4, 5},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3}
            };
  PllSplitList fst_splitlist = createSplitList(fst_part1s);

  std::vector<std::vector<size_t>> snd_part1s = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList snd_splitlist = createSplitList(snd_part1s);

  std::vector<std::vector<size_t>> trd_part1s = {
                { 0, 2, 3, 4, 5, 6, 7},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList trd_splitlist = createSplitList(trd_part1s);

  std::vector<std::vector<size_t>> delta12 = {
                { 0, 1, 2, 3, 4, 5},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3},
                { 0, 2, 3, 6}
            };
  PllSplitList delta12_splitlist = createSplitList(delta12);

  std::vector<std::vector<size_t>> delta13 = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 1, 2, 3, 4, 5},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3},
                { 0, 2, 3, 4, 5, 6, 7},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList delta13_splitlist = createSplitList(delta13);

  std::vector<std::vector<size_t>> delta23 = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6, 7},
  };
  PllSplitList delta23_splitlist = createSplitList(delta23);

  ASSERT_EQ(fst_splitlist.rfDistance(snd_splitlist), 6);
  ASSERT_EQ(snd_splitlist.rfDistance(fst_splitlist), 6);

  ASSERT_EQ(fst_splitlist.rfDistance(trd_splitlist), 9);
  ASSERT_EQ(trd_splitlist.rfDistance(fst_splitlist), 9);

  ASSERT_EQ(snd_splitlist.rfDistance(trd_splitlist), 3);
  ASSERT_EQ(trd_splitlist.rfDistance(snd_splitlist), 3);

  TestUtil::split_lists_eq(fst_splitlist.symmetricDifference(snd_splitlist), delta12_splitlist);
  TestUtil::split_lists_eq(snd_splitlist.symmetricDifference(fst_splitlist), delta12_splitlist);

  TestUtil::split_lists_eq(fst_splitlist.symmetricDifference(trd_splitlist), delta13_splitlist);
  TestUtil::split_lists_eq(trd_splitlist.symmetricDifference(fst_splitlist), delta13_splitlist);

  TestUtil::split_lists_eq(snd_splitlist.symmetricDifference(trd_splitlist), delta23_splitlist);
  TestUtil::split_lists_eq(trd_splitlist.symmetricDifference(snd_splitlist), delta23_splitlist);


}
