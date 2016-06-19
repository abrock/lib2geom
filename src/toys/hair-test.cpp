#include <gtest/gtest.h>

#include <toys/hair.h>

TEST(mod_functions, modu_dist) {
    Path a, b;
    Hair h(a,b);
    EXPECT_EQ(h.moduloDist(2,2,10), 0u);
    EXPECT_EQ(h.moduloDist(2,3,10), 1u);
    EXPECT_EQ(h.moduloDist(3,2,10), 1u);
    EXPECT_EQ(h.moduloDist(0,9,10), 1u);
    EXPECT_EQ(h.moduloDist(9,0,10), 1u);
    EXPECT_EQ(h.moduloDist(1,9,10), 2u);
    EXPECT_EQ(h.moduloDist(9,1,10), 2u);
    EXPECT_EQ(h.moduloDist(1,8,10), 3u);
    EXPECT_EQ(h.moduloDist(8,1,10), 3u);
}

TEST(mod_functions, modu_direction) {
    Path a, b;
    Hair h(a,b);
    EXPECT_EQ( 1l, h.moduloDistDirection(2,3,10));
    EXPECT_EQ(-1l, h.moduloDistDirection(3,2,10));
    EXPECT_EQ( 1l, h.moduloDistDirection(9,0,10));
    EXPECT_EQ(-1l, h.moduloDistDirection(0,9,10));
    EXPECT_EQ( 1l, h.moduloDistDirection(2,5,10));
    EXPECT_EQ(-1l, h.moduloDistDirection(5,2,10));
    EXPECT_EQ( 1l, h.moduloDistDirection(9,2,10));
    EXPECT_EQ(-1l, h.moduloDistDirection(2,9,10));
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}