#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

/**
 * Test is_symmetric() method.
 */
TEST(MatrixCheckerTest, is_symmetric_test)
{
    Matrix A(10, 10);
    // zero matrix should be regarded as symmetric.
    EXPECT_TRUE(A.is_symmetric());
    EXPECT_TRUE(A.is_symmetric(1e-16));
    EXPECT_TRUE(A.is_symmetric(1e-300));

    A.randomize(0, 1);
    EXPECT_FALSE(A.is_symmetric()); // random matrix usually is not symmetric.
    EXPECT_TRUE(A.is_symmetric(2)); // lower threshold to be symmetric.
    A.to_symmetric("L");
    EXPECT_TRUE(A.is_symmetric());  // real symmetric matrix.
    EXPECT_TRUE(A.is_symmetric(-1e-16)); // happended to set it to negative threshold.

    // if customized threshold working?
    double thred = 1e-8;
    A(2,3) = A(3, 2) - thred;
    EXPECT_FALSE(A.is_symmetric(thred/10));
    EXPECT_TRUE(A.is_symmetric(thred * 10));

    // non-square matrix.
    Matrix B(2, 3);
    EXPECT_FALSE(B.is_symmetric());
    EXPECT_FALSE(B.is_symmetric(100));
    EXPECT_FALSE(B.is_symmetric(1e-300));
}

/**
 * Test is_diagonal() method.
 */
TEST(MatrixCheckerTest, is_diagonal_test)
{
    Matrix A(10, 10);
    // zero matrix should be regarded as diagonal.
    EXPECT_TRUE(A.is_diagonal());
    EXPECT_TRUE(A.is_diagonal(1e-16));
    EXPECT_TRUE(A.is_diagonal(1e-300));

    // random matrix usually is not diagonal.
    A.randomize(0, 1);
    EXPECT_FALSE(A.is_diagonal());
    EXPECT_TRUE(A.is_diagonal(2));  // lower threshold to be true.

    // first make it diagonal matrix.
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < i; j++) {
            A(i, j) = A(j, i) = 0.0;
        }
    }
    EXPECT_TRUE(A.is_diagonal());
    EXPECT_TRUE(A.is_diagonal(1e-300));
    EXPECT_TRUE(A.is_diagonal(-1e-300)); // happended to set thre to negative.

    // customized threshold
    double thred = 1e-8;
    A(1,2) += thred;
    EXPECT_FALSE(A.is_diagonal(thred/10));
    EXPECT_TRUE(A.is_diagonal(thred * 10));

    // non-square matrix.
    Matrix B(2, 3);
    EXPECT_FALSE(B.is_diagonal());
    EXPECT_FALSE(B.is_diagonal(1e-300));
}

/**
 * Test is_identity() method.
 */
TEST(MatrixCheckerTest, is_identity_test)
{
    Matrix A(10, 10);
    // zero matrix should be not regarded as diagonal.
    EXPECT_FALSE(A.is_identity());
    EXPECT_FALSE(A.is_identity(1e-16));
    EXPECT_FALSE(A.is_identity(1e-300));

    // random matrix usually is not identity.
    A.randomize(0, 1);
    EXPECT_FALSE(A.is_identity());
    //EXPECT_TRUE(A.is_identity(2));  // lower threshold to be true.

    // first make it diagonal matrix.
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < i; j++) {
            A(i, j) = A(j, i) = 0.0;
        }
    }
    EXPECT_FALSE(A.is_identity());
    EXPECT_FALSE(A.is_identity(1e-300));

    // real identity matrix.
    for (size_t i = 0; i < A.row(); i++) {
        A(i, i) = 1.0;
    }
    EXPECT_TRUE(A.is_identity());
    EXPECT_TRUE(A.is_identity(1e-300));
    EXPECT_TRUE(A.is_identity(-1e-300)); // happended to set thred to negative.

    // customized threshold
    double thred = 1e-8;
    A(1,1) += thred;
    EXPECT_FALSE(A.is_identity(thred/10));
    EXPECT_TRUE(A.is_identity(thred * 10));

    // non-square matrix.
    Matrix B(2, 3);
    EXPECT_FALSE(B.is_identity());
    EXPECT_FALSE(B.is_identity(1e-300));
}

/**
 * Test is_zeros() method.
 */
TEST(MatrixCheckerTest, is_zeros_test)
{
    Matrix A(10, 8);
    // real zero matrix
    EXPECT_TRUE(A.is_zeros());
    EXPECT_TRUE(A.is_zeros(1e-16));
    EXPECT_TRUE(A.is_zeros(1e-300));

    // customized threshold
    double thred = 1e-8;
    A(1,1) += thred;
    EXPECT_FALSE(A.is_zeros(thred/10));
    EXPECT_TRUE(A.is_zeros(thred * 10));

    // random matrix usually is not zero.
    A.randomize(0, 1);
    EXPECT_FALSE(A.is_zeros());
    //EXPECT_TRUE(A.is_zeros(2));  // lower threshold to be true.
}

/**
 * Test is_equal_to() method.
 */
TEST(MatrixCheckerTest, is_equal_to_test)
{
    Matrix A(10, 8);
    Matrix B(10, 8);
    // real zero matrix should be equal.
    EXPECT_TRUE(A.is_equal_to(B));
    EXPECT_TRUE(A.is_equal_to(B, 1e-16));
    EXPECT_TRUE(A.is_equal_to(B, 1e-300));

    // customized threshold
    double thred = 1e-8;
    A(1,1) += thred;
    EXPECT_FALSE(A.is_equal_to(B, thred/10));
    EXPECT_TRUE(A.is_equal_to(B, thred * 10));

    // random matrix usually is not equal to each other.
    A.randomize(0, 1);
    B.randomize(0, 1);
    EXPECT_FALSE(A.is_equal_to(B));
}

TEST(MatrixCheckerTest, is_data_stored_outside)
{
    vector<double> data = {1, 2, 3, 4};
    // special case, empty matrix should reported with false.
    Matrix A;
    EXPECT_FALSE(A.is_data_stored_outside());
    Matrix B(0, 0);
    EXPECT_FALSE(B.is_data_stored_outside());

    Matrix C(2, 2, data, matrix::Matrix::CopyType::kShallowCopy);
    EXPECT_TRUE(C.is_data_stored_outside());
    Matrix D(2, 2, data, matrix::Matrix::CopyType::kDeepCopy);
    EXPECT_FALSE(D.is_data_stored_outside());
    Matrix E(2, 2);
    EXPECT_FALSE(E.is_data_stored_outside());
}
