#include "gtest/gtest.h"
#include "../../../src/PllSplits.hpp"
#include "../../../src/io/FileReader.hpp"
class PllSplitTest : public testing::Test {
protected:
  PllSplit createSplit(size_t split_len, std::vector<size_t> part1) {
    if (part1[0] != 0) throw "In every split, 0 must be in Partition 1, hence it must hold that part1[0]==0";
    auto split_bits = (pll_split_t)calloc(split_len, sizeof(pll_split_base_t));
    size_t major_idx;
    size_t minor_idx;
    for(size_t tip : part1){
      major_idx = tip / (sizeof(pll_split_base_t) * 8);
      minor_idx = tip % (sizeof(pll_split_base_t) * 8);
      split_bits[major_idx] |= (1 << minor_idx);
    }
    return PllSplit(split_bits, split_len);
  }

  PllSplitList createSplitList(size_t split_len, std::vector<std::vector<size_t>> part1s){
    std::vector<PllSplit> splits;
    for (auto part1 : part1s){
      splits.emplace_back(createSplit(split_len, part1));
    }
    std::sort(splits.begin(), splits.end());
    return PllSplitList(splits);
  }

  void  split_lists_eq(PllSplitList l1, PllSplitList l2){
    std::vector<PllSplit> splits1 = l1. getSplits();
    std::vector<PllSplit> splits2 = l2.getSplits();
    EXPECT_EQ(splits1.size(), splits2.size());
    for(size_t i = 0; i < splits1.size(); ++i){
      EXPECT_EQ(splits1[i], splits2[i]);
    }

  }

};


TEST_F(PllSplitTest, test_bitExtract) {
    // tips in the partition 1
    //0 must never be in part1 due to normalization
    std::vector<size_t> part1 {0, 2, 4, 8, 19, 45, 63};
    size_t split_len = 2;
    PllSplit split = createSplit(split_len, part1);
    std::vector<bool> split_representation(split_len*32);
    for(size_t tip : part1) {
      split_representation[tip] = true;
    }
    for(size_t i = 0; i < split_len*32; i++){
      EXPECT_EQ(split_representation[i], split.bitExtract(i));
    }
}


TEST_F(PllSplitTest, test_popcount) {
  // tips in the partition 1 (0 is always in, does not need to be added)
  std::vector<size_t> part1 {0, 2, 4, 8, 19, 45, 63};
  size_t split_len = 2;
  PllSplit split = createSplit(split_len, part1);
  EXPECT_EQ(part1.size(), split.popcount());
}


TEST_F(PllSplitTest, test_operators) {
  size_t split_len = 2;
  // tips in the partition 1
  // contains 0 impicitly
  std::vector<size_t> fst_part1 {0, 2, 4, 8, 19, 45, 63};
  PllSplit fst_split = createSplit(split_len, fst_part1);
  std::vector<size_t> snd_part1 {0, 1, 4, 8, 19, 45, 63};
  PllSplit snd_split = createSplit(split_len, snd_part1);
  ASSERT_TRUE(snd_split < fst_split);

  std::vector<size_t> trd_part1 {0, 2, 4, 8, 19, 45, 63};
  PllSplit trd_split = createSplit(split_len, trd_part1);
  ASSERT_TRUE(trd_split == fst_split);
  ASSERT_TRUE(snd_split < fst_split);
}

TEST_F(PllSplitTest, test_list_constructor) {
  size_t split_len = 2;
  std::vector<PllSplit> splits;
  std::vector<size_t> fst_part1 {0, 2, 4, 8, 19, 45, 63};
  splits.emplace_back(createSplit(split_len, fst_part1));
  std::vector<size_t> snd_part1 {0, 1, 4, 8, 19, 45, 63};
  splits.emplace_back(createSplit(split_len, snd_part1));
  std::vector<size_t> trd_part1 {0, 19, 29, 39};
  splits.emplace_back(createSplit(split_len, trd_part1));
  PllSplitList split_list = PllSplitList(splits);
  split_lists_eq(splits, split_list.getSplits());

}

TEST_F(PllSplitTest, test_tree_constructor) {
  PllSplitList split_list_from_tree = PllSplitList(PllTree("((a1, a2), (b1,b2), (c, (d1, d2)));"));
  size_t split_len = 1;
  std::vector<std::vector<size_t>> expected {
                { 0, 1},
                { 0, 1, 2, 3},
                { 0, 1, 4, 5, 6},
                { 0, 1, 2, 3, 4},
  };
  PllSplitList expected_splitlist = createSplitList(split_len, expected);
  split_lists_eq(expected_splitlist, split_list_from_tree);

}


TEST_F(PllSplitTest, test_difference) {
  size_t split_len = 1;
  // tips in the partition 1
  // contains 0 impicitly
  std::vector<std::vector<size_t>> fst_part1s {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 1, 2, 3, 4, 5},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3}
            };
  PllSplitList fst_splitlist = createSplitList(split_len, fst_part1s);

  std::vector<std::vector<size_t>> snd_part1s {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList snd_splitlist = createSplitList(split_len, snd_part1s);

  std::vector<std::vector<size_t>> trd_part1s {
                { 0, 2, 3, 4, 5, 6, 7},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 2, 3, 6}
            };
  PllSplitList trd_splitlist = createSplitList(split_len, trd_part1s);

  std::vector<std::vector<size_t>> delta12 {
                { 0, 1, 2, 3, 4, 5},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3},
                { 0, 2, 3, 6}
            };
  PllSplitList delta12_splitlist = createSplitList(split_len, delta12);

  std::vector<std::vector<size_t>> delta13 {
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
  PllSplitList delta13_splitlist = createSplitList(split_len, delta13);

  std::vector<std::vector<size_t>> delta23 {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6, 7},
  };
  PllSplitList delta23_splitlist = createSplitList(split_len, delta23);

  ASSERT_EQ(fst_splitlist.rfDistance(snd_splitlist), 6);
  ASSERT_EQ(snd_splitlist.rfDistance(fst_splitlist), 6);

  ASSERT_EQ(fst_splitlist.rfDistance(trd_splitlist), 9);
  ASSERT_EQ(trd_splitlist.rfDistance(fst_splitlist), 9);

  ASSERT_EQ(snd_splitlist.rfDistance(trd_splitlist), 3);
  ASSERT_EQ(trd_splitlist.rfDistance(snd_splitlist), 3);

  split_lists_eq(fst_splitlist.symmetricDifference(snd_splitlist), delta12_splitlist);
  split_lists_eq(snd_splitlist.symmetricDifference(fst_splitlist), delta12_splitlist);

  split_lists_eq(fst_splitlist.symmetricDifference(trd_splitlist), delta13_splitlist);
  split_lists_eq(trd_splitlist.symmetricDifference(fst_splitlist), delta13_splitlist);

  split_lists_eq(snd_splitlist.symmetricDifference(trd_splitlist), delta23_splitlist);
  split_lists_eq(trd_splitlist.symmetricDifference(snd_splitlist), delta23_splitlist);


}
