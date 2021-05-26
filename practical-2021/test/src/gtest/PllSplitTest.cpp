#include "gtest/gtest.h"
#include "../../../src/datastructures/PllSplits.hpp"
#include "../../../src/io/TreeReader.hpp"
#include "../TestUtil.hpp"


class PllSplitTest : public testing::Test {

protected:

  void split_vector_eq(const std::vector<PllSplit>& l1, const std::vector<PllSplit>& l2) {
    EXPECT_EQ(l1.size(), l2.size());
    for(size_t i = 0; i < l1.size(); ++i){
      EXPECT_EQ(l1[i], l2[i]);
    }
  }
  void constructor_eq() {
    PllTree test_tree = PllTree("((a1, a2), (b1,b2), (c, (d1, d2)));");
    PllSplit::setTipCount(test_tree.getTipCount());
    PllSplitList split_list_from_tree = PllSplitList(test_tree);
    std::vector<std::vector<size_t>> expected = {
                { 0, 1},
                { 0, 1, 2, 3},
                { 0, 1, 4, 5, 6},
                { 0, 1, 2, 3, 4},
    };
    PllSplitList expected_splitlist = TestUtil::createSplitList(expected);
    ASSERT_EQ(expected_splitlist, split_list_from_tree);
  }
};


TEST_F(PllSplitTest, test_bitExtract) {
    PllSplit::setTipCount(64);
    std::vector<size_t> part1 = {0, 2, 4, 8, 19, 45, 63};
    PllSplit split = TestUtil::TestUtil::createSplit(part1);
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
  PllSplit::setTipCount(64);
  // tips in the partition 1 (0 needs to be in)
  std::vector<size_t> part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit split = TestUtil::createSplit(part1);
  EXPECT_EQ(part1.size(), split.popcount());
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

TEST_F(PllSplitTest, test_list_constructor) {
  PllSplit::setTipCount(64);
  std::vector<PllSplit> splits;
  std::vector<size_t> fst_part1 = {0, 2, 4, 8, 19, 45, 63};
  splits.emplace_back(TestUtil::createSplit(fst_part1));
  std::vector<size_t> snd_part1 = {0, 1, 4, 8, 19, 45, 63};
  splits.emplace_back(TestUtil::createSplit(snd_part1));
  std::vector<size_t> trd_part1 = {0, 19, 29, 39};
  splits.emplace_back(TestUtil::createSplit(trd_part1));
  PllSplitList split_list = PllSplitList(splits);
  split_vector_eq(splits, split_list.getSplits());
  for(PllSplit split : splits){
    free(split());
  }

}

TEST_F(PllSplitTest, test_tree_constructor) {
  constructor_eq();

}
TEST_F(PllSplitTest, test_tree_constructor_iterative) {
  for(size_t i = 0; i < 100; ++i) {
    std::cout << "Running test_tree_constructor with iteration: " << i << "\n";
    constructor_eq();
  }

}


TEST_F(PllSplitTest, test_difference) {
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

  std::vector<std::vector<size_t>> delta12 = {
                { 0, 1, 2, 3, 4, 5},
                { 0, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 6},
                { 0, 1, 2, 3, 4},
                { 0, 1, 2, 3},
                { 0, 2, 3, 6}
            };
  PllSplitList delta12_splitlist = TestUtil::createSplitList(delta12);

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
  PllSplitList delta13_splitlist = TestUtil::createSplitList(delta13);

  std::vector<std::vector<size_t>> delta23 = {
                { 0, 1, 2, 3, 4, 5, 6, 7},
                { 0, 1, 2, 3, 4, 5, 6},
                { 0, 2, 3, 4, 5, 6, 7},
  };
  PllSplitList delta23_splitlist = TestUtil::createSplitList(delta23);

  ASSERT_EQ(fst_splitlist.rfDistance(snd_splitlist), 6);
  ASSERT_EQ(snd_splitlist.rfDistance(fst_splitlist), 6);

  ASSERT_EQ(fst_splitlist.rfDistance(trd_splitlist), 9);
  ASSERT_EQ(trd_splitlist.rfDistance(fst_splitlist), 9);

  ASSERT_EQ(snd_splitlist.rfDistance(trd_splitlist), 3);
  ASSERT_EQ(trd_splitlist.rfDistance(snd_splitlist), 3);

  ASSERT_EQ(fst_splitlist.symmetricDifference(snd_splitlist), delta12_splitlist);
  ASSERT_EQ(snd_splitlist.symmetricDifference(fst_splitlist), delta12_splitlist);

  ASSERT_EQ(fst_splitlist.symmetricDifference(trd_splitlist), delta13_splitlist);
  ASSERT_EQ(trd_splitlist.symmetricDifference(fst_splitlist), delta13_splitlist);

  ASSERT_EQ(snd_splitlist.symmetricDifference(trd_splitlist), delta23_splitlist);
  ASSERT_EQ(trd_splitlist.symmetricDifference(snd_splitlist), delta23_splitlist);

}



TEST_F(PllSplitTest, test_intersectcount) {
    PllSplit::setTipCount(10);
    std::vector<size_t> part1_a = {0, 2, 4, 9};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 1, 7, 9};
    PllSplit split_b = TestUtil::createSplit(part1_b);

    EXPECT_EQ(split_a.intersectionSize(split_a, 1, 1), 4);
    EXPECT_EQ(split_a.intersectionSize(split_a, 0, 0), 6);
    EXPECT_EQ(split_a.intersectionSize(split_a, 1, 0), 0);
    EXPECT_EQ(split_a.intersectionSize(split_a, 0, 1), 0);

    EXPECT_EQ(split_a.intersectionSize(split_b, 1, 1), 2);
    EXPECT_EQ(split_a.intersectionSize(split_b, 0, 0), 4);
    EXPECT_EQ(split_a.intersectionSize(split_b, 1, 0), 2);
    EXPECT_EQ(split_a.intersectionSize(split_b, 0, 1), 2);

    free(split_a());
    free(split_b());
}


TEST_F(PllSplitTest, test_unioncount) {
    PllSplit::setTipCount(6);
    std::vector<size_t> part1_a = {0, 1, 2, 3};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 1, 4, 5};
    PllSplit split_b = TestUtil::createSplit(part1_b);

    EXPECT_EQ(split_a.unionSize(split_a, 1, 1), 4);
    EXPECT_EQ(split_a.unionSize(split_a, 0, 0), 2);
    EXPECT_EQ(split_a.unionSize(split_a, 1, 0), 6);
    EXPECT_EQ(split_a.unionSize(split_a, 0, 1), 6);

    EXPECT_EQ(split_a.unionSize(split_b, 1, 1), 6);
    EXPECT_EQ(split_a.unionSize(split_b, 0, 0), 4);
    EXPECT_EQ(split_a.unionSize(split_b, 1, 0), 4);
    EXPECT_EQ(split_a.unionSize(split_b, 0, 1), 4);

    free(split_a());
    free(split_b());
}


TEST_F(PllSplitTest, test_compatible) {
    PllSplit::setTipCount(8);
    std::vector<size_t> part1_a = {0, 1, 2, 3};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 3};
    PllSplit split_b = TestUtil::createSplit(part1_b);
    std::vector<size_t> part1_c = {0, 3, 4};
    PllSplit split_c = TestUtil::createSplit(part1_c);
    EXPECT_TRUE(split_a.compatible(split_b));
    EXPECT_TRUE(split_b.compatible(split_a));
    EXPECT_FALSE(split_a.compatible(split_c));
    EXPECT_FALSE(split_c.compatible(split_a));
    EXPECT_TRUE(split_c.compatible(split_b));
    EXPECT_TRUE(split_b.compatible(split_c));

    free(split_a());
    free(split_b());
    free(split_c());
}


TEST_F(PllSplitTest, test_containsassubset) {
    PllSplit::setTipCount(6);
    std::vector<size_t> part1_a = {0, 1, 2, 3};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 1};
    PllSplit split_b = TestUtil::createSplit(part1_b);
    std::vector<size_t> part1_c = {0, 1, 4, 5};
    PllSplit split_c = TestUtil::createSplit(part1_c);

    EXPECT_TRUE(split_a.containsAsSubset(split_b, 1, 1));
    EXPECT_TRUE(split_b.containsAsSubset(split_a, 0, 0));
    EXPECT_TRUE(split_c.containsAsSubset(split_a, 1, 0));
    EXPECT_FALSE(split_a.containsAsSubset(split_c, 1, 1));

    free(split_a());
    free(split_b());
    free(split_c());
}


TEST_F(PllSplitTest, test_containsassubset2) {
    PllSplit::setTipCount(24);
    std::vector<size_t> part1_a = {0, 1, 5, 17, 22};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 1, 22};
    PllSplit split_b = TestUtil::createSplit(part1_b);

    EXPECT_TRUE(split_a.containsAsSubset(split_b, 1, 1));
    EXPECT_TRUE(split_b.containsAsSubset(split_a, 0, 0));
    EXPECT_FALSE(split_b.containsAsSubset(split_a, 1, 1));

    free(split_a());
    free(split_b());
}
