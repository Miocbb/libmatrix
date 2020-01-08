#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>
#include "utils.h"

using matrix::Matrix;
using std::vector;
using Eigen::MatrixXd;

/**
 * test functions to invert matrix.
 */
struct InvertTest: public ::testing::Test {
    Matrix A22_gen_mat;
    Matrix A22_gen_inv_mat;

    Matrix A22_sym_mat;
    Matrix A22_sym_inv_mat;

    virtual void SetUp() override
    {
        A22_gen_mat = Matrix(2, 2).randomize(0, 1);
        MatrixXd A22_gen_mxd = Matrix_to_MatrixXd(A22_gen_mat);
        MatrixXd A22_gen_inv_mxd = A22_gen_mxd.inverse();
        A22_gen_inv_mat = MatrixXd_to_Matrix(A22_gen_inv_mxd);

        A22_sym_mat = Matrix(2, 2).randomize(0, 1).to_symmetric("L");
        MatrixXd A22_sym_mxd = Matrix_to_MatrixXd(A22_sym_mat);
        MatrixXd A22_sym_inv_mxd = A22_sym_mxd.inverse();
        A22_sym_inv_mat = MatrixXd_to_Matrix(A22_sym_inv_mxd);
    }

    virtual void TearDown() override {}
};

/**
 * testing function using dgetri subroutine.
 */
TEST_F(InvertTest, general_matrix_dgetri_test)
{
    Matrix A;

    // invert a general matrix
    A = A22_gen_mat;
    matrix::invert_gen_matrix_dgetri(A);
    EXPECT_TRUE(A.is_equal_to(A22_gen_inv_mat, 1e-10));

    // invert a sym matrix
    A = A22_sym_mat;
    matrix::invert_gen_matrix_dgetri(A);
    EXPECT_TRUE(A.is_equal_to(A22_sym_inv_mat, 1e-10));
}

///**
// * testing function using dpotri subroutine.
// */
//TEST_F(InvertTest, spd_matrix_dpotri_test)
//{
//    Matrix A;
//
//    // A22_sym_mat is a symmetric random matrix.
//    // most of case, it should be indefinite, so it should give false results.
//    A = A22_sym_mat;
//    matrix::invert_spd_matrix_dpotri("U", A);
//    EXPECT_FALSE(A.is_equal_to(A22_sym_inv_mat));
//    A = A22_sym_mat;
//    matrix::invert_spd_matrix_dpotri("L", A);
//    EXPECT_FALSE(A.is_equal_to(A22_sym_inv_mat));
//}

/**
 * testing function using dsytri subroutine.
 */
TEST_F(InvertTest, symmetric_matrix_dsytri_test)
{
    Matrix A;

    // invert a sym matrix
    A = A22_sym_mat;
    A(1, 0) = 999;
    matrix::invert_sym_matrix_dsytri("U", A);
    EXPECT_TRUE(A.is_equal_to(A22_sym_inv_mat, 1e-10));
    A = A22_sym_mat;
    A(0, 1) = 999;
    matrix::invert_sym_matrix_dsytri("L", A);
    EXPECT_TRUE(A.is_equal_to(A22_sym_inv_mat, 1e-10));
}

/**
 * testing function using dsytri_rook subroutine.
 */
TEST_F(InvertTest, symmetric_matrix_dsytri_rook_test)
{
    Matrix A;

    // invert a sym matrix
    A = A22_sym_mat;
    A(1, 0) = 999;
    matrix::invert_sym_matrix_dsytri_rook("U", A);
    EXPECT_TRUE(A.is_equal_to(A22_sym_inv_mat, 1e-10));
    A = A22_sym_mat;
    A(0, 1) = 999;
    matrix::invert_sym_matrix_dsytri_rook("L", A);
    EXPECT_TRUE(A.is_equal_to(A22_sym_inv_mat, 1e-10));
}
