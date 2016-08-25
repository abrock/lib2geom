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

TEST(path_functions, curveLength) {
    Point a(0,0);
    Point b(1,0);
    Point c(0,1);
    BezierCurve _straightLine = BezierCurveN<1>(a,b);
    Path straightLine;
    straightLine.append(_straightLine);
    PathTime start, stop;
    stop += 1;
    EXPECT_NEAR(Hair::curveLength(straightLine, start, stop), 1.0, 1e-10);
    stop += -.5;
    EXPECT_NEAR(Hair::curveLength(straightLine, start, stop), 0.5, 1e-10);
    Path curve = straightLine;
    _straightLine = BezierCurveN<1>(b,c);
    curve.append(_straightLine);
    stop = start;
    stop += 2;
    EXPECT_NEAR(Hair::curveLength(curve, start, stop), 1.0 + M_SQRT2, 1e-10);
    stop += -.5;
    EXPECT_NEAR(Hair::curveLength(curve, start, stop), 1.0 + 0.5 * M_SQRT2, 1e-10);
    start += .5;
    EXPECT_NEAR(Hair::curveLength(curve, start, stop), 0.5 + 0.5 * M_SQRT2, 1e-10);
    stop += .25;
    EXPECT_NEAR(Hair::curveLength(curve, start, stop), 0.5 + 0.75 * M_SQRT2, 1e-10);
    start += .25;
    EXPECT_NEAR(Hair::curveLength(curve, start, stop), 0.25 + 0.75 * M_SQRT2, 1e-10);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    std::cout << "RUN_ALL_TESTS return value: " << RUN_ALL_TESTS() << std::endl;
}
