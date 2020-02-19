#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>


TEST(MatrixResizeTest, data_stored_inside)
{
    matrix::Matrix A;
    A.resize(2, 2);
    EXPECT_TRUE(A.is_zeros());

    A = {1, 2,
         3, 4};
    A.resize(1, 2);
    matrix::Matrix A_ref_1(1, 2);
    A_ref_1 = {1, 2};
    EXPECT_TRUE(A.is_equal_to(A_ref_1));

    A.resize(2, 2);
    matrix::Matrix A_ref_2(2, 2);
    A_ref_2 = {1, 2,
               0, 0};
    EXPECT_TRUE(A.is_equal_to(A_ref_2));
}

TEST(MatrixResizeTest, data_stored_outside)
{

    std::vector<double> d = {1, 2, 3, 4};
    matrix::Matrix A(2, 2, d, matrix::Matrix::kShallowCopy);
    A.resize(1, 2);
    matrix::Matrix A_ref_1(1, 2);
    A_ref_1 = {1, 2};
    EXPECT_TRUE(A.is_equal_to(A_ref_1));
    EXPECT_FALSE(A.is_data_stored_outside());

    d = {1, 2};
    matrix::Matrix A2(1, 2, d, matrix::Matrix::kShallowCopy);
    A2.resize(2, 2);
    matrix::Matrix A_ref_2(2, 2);
    A_ref_2 = {1, 2,
               0, 0};
    EXPECT_TRUE(A2.is_equal_to(A_ref_2));
    EXPECT_FALSE(A2.is_data_stored_outside());
}
