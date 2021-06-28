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

/*TEST_F(PllSplitTest, test_list_constructor) {
  PllSplit::setTipCount(64);
  std::vector<PllSplit*> splits;
  std::vector<size_t> fst_part1 = {0, 2, 4, 8, 19, 45, 63};
  PllSplit t1 = TestUtil::createSplit(fst_part1);



  splits.emplace_back(&t1);
  std::vector<size_t> snd_part1 = {0, 1, 4, 8, 19, 45, 63};
    PllSplit t2 = TestUtil::createSplit(snd_part1);
  splits.emplace_back(&t2);
  std::vector<size_t> trd_part1 = {0, 19, 29, 39};
    PllSplit t3 = TestUtil::createSplit(trd_part1);
  splits.emplace_back(&t3);
  PllSplitList split_list = PllSplitList(splits);
  split_vector_eq(splits, split_list.getSplits());
  for(PllSplit* split : splits){
    //free(split);
  }

}*/

/*TEST_F(PllSplitTest, test_tree_constructor) {
  constructor_eq();

}
TEST_F(PllSplitTest, test_tree_constructor_iterative) {
  for(size_t i = 0; i < 100; ++i) {
    //std::cout << "Running test_tree_constructor with iteration: " << i << "\n";
    constructor_eq();
  }

}*/


TEST_F(PllSplitTest, test_intersectcount) {
    PllSplit::setTipCount(10);
    std::vector<size_t> part1_a = {0, 2, 4, 9};
    PllSplit split_a = TestUtil::createSplit(part1_a);
    std::vector<size_t> part1_b = {0, 1, 7, 9};
    PllSplit split_b = TestUtil::createSplit(part1_b);

    EXPECT_EQ(split_a.intersectionSize(split_a, Block_A, Block_A), 4);
    EXPECT_EQ(split_a.intersectionSize(split_a, Block_B, Block_B), 6);
    EXPECT_EQ(split_a.intersectionSize(split_a, Block_A, Block_B), 0);
    EXPECT_EQ(split_a.intersectionSize(split_a, Block_B, Block_A), 0);

    EXPECT_EQ(split_a.intersectionSize(split_b, Block_A, Block_A), 2);
    EXPECT_EQ(split_a.intersectionSize(split_b, Block_B, Block_B), 4);
    EXPECT_EQ(split_a.intersectionSize(split_b, Block_A, Block_B), 2);
    EXPECT_EQ(split_a.intersectionSize(split_b, Block_B, Block_A), 2);

    free(split_a());
    free(split_b());
}





/*TEST_F(PllSplitTest, test_compatible) {
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
//This test is based on Luises 8 taxa tree which brings compatibility to its limit
//Since all the splits have been extracted from one tree they should all be compatible
TEST_F(PllSplitTest, test_compatible_8taxa) {
  PllSplit::setTipCount(8);
  PllSplit s1 = TestUtil::createSplit({0, 1}); // AB
  PllSplit s2 = TestUtil::createSplit({0, 1, 4, 5, 6, 7}); //ABEFGH
  PllSplit s3 = TestUtil::createSplit({0, 1, 2, 3, 6, 7}); //ABCDGH
  PllSplit s4 = TestUtil::createSplit({0, 1, 2, 3, 4, 5}); //ABCDEF
  PllSplit s5 = TestUtil::createSplit({0, 1, 2, 3});  //ABCD
  EXPECT_TRUE(s1.compatible(s2));
  EXPECT_TRUE(s1.compatible(s3));
  EXPECT_TRUE(s1.compatible(s4));
  EXPECT_TRUE(s1.compatible(s5));
  EXPECT_TRUE(s2.compatible(s3));
  EXPECT_TRUE(s2.compatible(s4));
  EXPECT_TRUE(s2.compatible(s5));
  EXPECT_TRUE(s3.compatible(s4));
  EXPECT_TRUE(s3.compatible(s5));
  EXPECT_TRUE(s4.compatible(s5));

  free(s1());
  free(s2());
  free(s3());
  free(s4());
  free(s5());

}*/





//These test cases were added because of the deprecated method of bitmask
TEST_F(PllSplitTest, popcount_2registers_nonfull) {
  PllSplit::setTipCount(36);
  PllSplit test_split = TestUtil::createSplit({0, 31, 33, 35});
  EXPECT_EQ(test_split.popcount(), 4);
  free(test_split());
}
TEST_F(PllSplitTest, intersectionsize_2registers_nonfull) {
  PllSplit::setTipCount(36);
  PllSplit test_split = TestUtil::createSplit({0, 31, 33, 35});
  EXPECT_EQ(test_split.intersectionSize(test_split, Block_A, Block_A), 4);
  EXPECT_EQ(test_split.intersectionSize(test_split, Block_B, Block_A), 0);
  EXPECT_EQ(test_split.intersectionSize(test_split, Block_A, Block_B), 0);
  EXPECT_EQ(test_split.intersectionSize(test_split, Block_B, Block_B), 32);
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
