#include "gtest/gtest.h"

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(testMath, myCubeTest)
{
    EXPECT_EQ(1, 1);
    std::cout << " HELLO WORLD" << std::endl;
}
