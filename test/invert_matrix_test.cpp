#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>
#include <matrix/matrix.h>
#include <random>
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
        A22_gen_mat.resize(2, 2);
        A22_gen_mat.randomize(0, 1);
        MatrixXd A22_gen_mxd = Matrix_to_MatrixXd(A22_gen_mat);
        MatrixXd A22_gen_inv_mxd = A22_gen_mxd.inverse();
        A22_gen_inv_mat = MatrixXd_to_Matrix(A22_gen_inv_mxd);

        A22_sym_mat.resize(2, 2);
        A22_sym_mat.randomize(0, 1);
        A22_sym_mat.to_symmetric("L");
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


TEST_F(InvertTest, spd_matrix_dpotri_test)
{
    Matrix Q(3, 3);
    matrix::set_matrix_random_orthogonal(Q);
    Matrix DiagM(3, 3);
    Matrix DiagM_inv(3, 3);

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<> dist6(1, 10);
    std::vector<double> rand_vec;
    for (size_t i = 0; i < Q.row(); ++i)
        rand_vec.push_back(dist6(rng));
    std::sort(rand_vec.begin(), rand_vec.end());
    for (size_t i = 0; i < Q.row(); ++i) {
        DiagM(i, i) = rand_vec[i];
        DiagM_inv(i, i) = 1.0/rand_vec[i];
    }

    // make spd matrix.
    Matrix A(3, 3);
    Matrix A_inv_ref(3, 3);
    std::vector<double> eigV(3);
    matrix::mult_dgemm_ATBA(Q, DiagM, A);
    matrix::mult_dgemm_ATBA(Q, DiagM_inv, A_inv_ref);

    // do diagonalization with Upper part.
    Matrix A_inv_calc = A;
    matrix::invert_spd_matrix_dpotri("U", A_inv_calc);
    EXPECT_TRUE(A_inv_calc.is_equal_to(A_inv_ref));

    // do diagonalization with Lower part.
    A_inv_calc = A;
    matrix::invert_spd_matrix_dpotri("L", A_inv_calc);
    EXPECT_TRUE(A_inv_calc.is_equal_to(A_inv_ref));
}
