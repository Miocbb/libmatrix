#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <Eigen/Dense>
#include <vector>
#include "utils.h"

using matrix::Matrix;
using matrix::mult_dgemm;
using Eigen::MatrixXd;

struct DgemmTest: public ::testing::Test {
    Matrix A22_gen_m;
    Matrix A23_gen_m;
    Matrix A32_gen_m;
    Matrix A33_gen_m;
    Matrix A100_m;
    MatrixXd A22_gen_mxd;
    MatrixXd A23_gen_mxd;
    MatrixXd A32_gen_mxd;
    MatrixXd A33_gen_mxd;
    MatrixXd A100_mxd;

    Matrix C22_m;
    Matrix C23_m;
    Matrix C32_m;
    Matrix C33_m;
    Matrix C100_m;

    MatrixXd C22_mxd;
    MatrixXd C23_mxd;
    MatrixXd C32_mxd;
    MatrixXd C33_mxd;
    MatrixXd C100_mxd;

    virtual void SetUp() override
    {
        A22_gen_m = Matrix(2, 2).randomize(0, 1);
        A23_gen_m = Matrix(2, 3).randomize(0, 1);
        A32_gen_m = Matrix(3, 2).randomize(0, 1);
        A33_gen_m = Matrix(3, 3).randomize(0, 1);
        A100_m = Matrix(100, 100).randomize(0, 100);

        A22_gen_mxd = Matrix_to_mrixXd(A22_gen_m);
        A23_gen_mxd = Matrix_to_mrixXd(A23_gen_m);
        A32_gen_mxd = Matrix_to_mrixXd(A32_gen_m);
        A33_gen_mxd = Matrix_to_mrixXd(A33_gen_m);
        A100_mxd = Matrix_to_mrixXd(A100_m);

        C22_m = Matrix(2, 2);
        C23_m = Matrix(2, 3);
        C32_m = Matrix(3, 2);
        C33_m = Matrix(3, 3);
        C100_m = Matrix(100, 100);

        C22_mxd = Matrix_to_mrixXd(C22_m);
        C23_mxd = Matrix_to_mrixXd(C23_m);
        C32_mxd = Matrix_to_mrixXd(C32_m);
        C33_mxd = Matrix_to_mrixXd(C33_m);
        C100_mxd = Matrix_to_mrixXd(C100_m);
    }

    virtual void TearDown() override
    {

    }
};

TEST_F(DgemmTest, NxN_test)
{
    // A22 x A22
    mult_dgemm(1.0, A22_gen_m, "N",
               A22_gen_m, "N",
               0.0, C22_m);
    C22_mxd = A22_gen_mxd * A22_gen_mxd;
    check_data_equality_with_EigenMatrix(C22_mxd, C22_m);

    // A22 x A23
    mult_dgemm(1.0, A22_gen_m, "N",
               A23_gen_m, "N",
               0.0, C23_m);
    C23_mxd = A22_gen_mxd * A23_gen_mxd;
    check_data_equality_with_EigenMatrix(C23_mxd, C23_m);

    // A32 x A22
    mult_dgemm(1.0, A32_gen_m, "N",
               A22_gen_m, "N",
               0.0, C32_m);
    C32_mxd = A32_gen_mxd * A22_gen_mxd;
    check_data_equality_with_EigenMatrix(C32_mxd, C32_m);

    // A32 x A23
    mult_dgemm(1.0, A32_gen_m, "N",
               A23_gen_m, "N",
               0.0, C33_m);
    C33_mxd = A32_gen_mxd * A23_gen_mxd;
    check_data_equality_with_EigenMatrix(C33_mxd, C33_m);

    // A100 x A100
    mult_dgemm(1.0, A100_m, "N",
               A100_m, "N",
               0.0, C100_m);
    C100_mxd = A100_mxd * A100_mxd;
    check_data_equality_with_EigenMatrix(C100_mxd, C100_m);
}

TEST_F(DgemmTest, NxT_test)
{
    // A22 x A22.T
    mult_dgemm(1.0, A22_gen_m, "N",
               A22_gen_m, "T",
               0.0, C22_m);
    C22_mxd = A22_gen_mxd * A22_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C22_mxd, C22_m);

    // A22 x A32.T
    mult_dgemm(1.0, A22_gen_m, "N",
               A32_gen_m, "T",
               0.0, C23_m);
    C23_mxd = A22_gen_mxd * A32_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C23_mxd, C23_m);

    // A32 x A22.T
    mult_dgemm(1.0, A32_gen_m, "N",
               A22_gen_m, "T",
               0.0, C32_m);
    C32_mxd = A32_gen_mxd * A22_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C32_mxd, C32_m);

    // A32 x A32.T
    mult_dgemm(1.0, A32_gen_m, "N",
               A32_gen_m, "T",
               0.0, C33_m);
    C33_mxd = A32_gen_mxd * A32_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C33_mxd, C33_m);

    // A100 x A100.T
    mult_dgemm(1.0, A100_m, "N",
               A100_m, "T",
               0.0, C100_m);
    C100_mxd = A100_mxd * A100_mxd.transpose();
    check_data_equality_with_EigenMatrix(C100_mxd, C100_m);
}

TEST_F(DgemmTest, TxN_test)
{
    // A22.T x A22
    mult_dgemm(1.0, A22_gen_m, "T",
               A22_gen_m, "N",
               0.0, C22_m);
    C22_mxd = A22_gen_mxd.transpose() * A22_gen_mxd;
    check_data_equality_with_EigenMatrix(C22_mxd, C22_m);

    // A22.T x A23
    mult_dgemm(1.0, A22_gen_m, "T",
               A23_gen_m, "N",
               0.0, C23_m);
    C23_mxd = A22_gen_mxd.transpose() * A23_gen_mxd;
    check_data_equality_with_EigenMatrix(C23_mxd, C23_m);

    // A23.T x A22
    mult_dgemm(1.0, A23_gen_m, "T",
               A22_gen_m, "N",
               0.0, C32_m);
    C32_mxd = A23_gen_mxd.transpose() * A22_gen_mxd;
    check_data_equality_with_EigenMatrix(C32_mxd, C32_m);

    // A32.T x A32
    mult_dgemm(1.0, A32_gen_m, "T",
               A32_gen_m, "N",
               0.0, C22_m);
    C22_mxd = A32_gen_mxd.transpose() * A32_gen_mxd;
    check_data_equality_with_EigenMatrix(C22_mxd, C22_m);

    // A100.T x A100
    mult_dgemm(1.0, A100_m, "T",
               A100_m, "N",
               0.0, C100_m);
    C100_mxd = A100_mxd.transpose() * A100_mxd;
    check_data_equality_with_EigenMatrix(C100_mxd, C100_m);
}

TEST_F(DgemmTest, TxT_test)
{
    // A22.T x A22.T
    mult_dgemm(1.0, A22_gen_m, "T",
               A22_gen_m, "T",
               0.0, C22_m);
    C22_mxd = A22_gen_mxd.transpose() * A22_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C22_mxd, C22_m);

    // A22.T x A32.T
    mult_dgemm(1.0, A22_gen_m, "T",
               A32_gen_m, "T",
               0.0, C23_m);
    C23_mxd = A22_gen_mxd.transpose() * A32_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C23_mxd, C23_m);

    // A23.T x A22.T
    mult_dgemm(1.0, A23_gen_m, "T",
               A22_gen_m, "T",
               0.0, C32_m);
    C32_mxd = A23_gen_mxd.transpose() * A22_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C32_mxd, C32_m);

    // A32.T x A23.T
    mult_dgemm(1.0, A32_gen_m, "T",
               A23_gen_m, "T",
               0.0, C22_m);
    C22_mxd = A32_gen_mxd.transpose() * A23_gen_mxd.transpose();
    check_data_equality_with_EigenMatrix(C22_mxd, C22_m);

    // A100.T x A100.T
    mult_dgemm(1.0, A100_m, "T",
               A100_m, "T",
               0.0, C100_m);
    C100_mxd = A100_mxd.transpose() * A100_mxd.transpose();
    check_data_equality_with_EigenMatrix(C100_mxd, C100_m);
}
