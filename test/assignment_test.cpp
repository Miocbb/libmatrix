#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

TEST(CopyAssignmentTest, test)
{
    vector<double> data = {1, 2, 3, 4};
    Matrix A(2, 2, data, matrix::Matrix::CopyType::kDeepCopy);
    Matrix B;
    B = A;
    EXPECT_TRUE(B.is_equal_to(A));
    B(1, 1) = 999;
    EXPECT_FALSE(B.is_equal_to(A));

    Matrix C;
    C = B = A;
    EXPECT_TRUE(B.is_equal_to(A));
    EXPECT_TRUE(C.is_equal_to(A));
    EXPECT_TRUE(C.is_equal_to(B));
    B(1, 1) = 999;
    EXPECT_FALSE(B.is_equal_to(A));
    EXPECT_FALSE(B.is_equal_to(C));
    EXPECT_TRUE(C.is_equal_to(A));
    C(1, 1) = 111;
    EXPECT_FALSE(B.is_equal_to(A));
    EXPECT_FALSE(C.is_equal_to(A));
    EXPECT_FALSE(C.is_equal_to(B));

    C = C;
    EXPECT_TRUE(C.is_equal_to(C));
}
