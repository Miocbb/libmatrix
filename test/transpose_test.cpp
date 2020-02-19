#include <algorithm>
#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <random>
#include <vector>

using std::vector;
using matrix::Matrix;

TEST(TransposeTest, general_test)
{
    Matrix A22(2, 2);
    Matrix A22T(2, 2);
    A22 =   {1, 2,
             3, 4};
    A22T =  {1, 3,
             2, 4};
    A22.transpose();
    EXPECT_TRUE(A22.is_equal_to(A22T));

    Matrix A23(2, 3);
    Matrix A23T(3, 2);
    A23 =   {1, 2, 3,
             4, 5, 6};
    A23T =  {1, 4,
             2, 5,
             3, 6};
    A23.transpose();
    EXPECT_TRUE(A23.is_equal_to(A23T));
}
