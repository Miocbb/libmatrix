#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "utils.h"

using matrix::Matrix;
using Eigen::MatrixXd;
using matrix::mult_dscal;
using matrix::mult_dscal_to;

/**
 * Test matrix scale method.
 */
TEST(MatrixScaleTest, method_scale_test)
{
    Matrix A(10, 100);
    Matrix B(10, 100);
    Matrix A_copy;
    Matrix A_copy2;
    A_copy2 = A_copy = A;
    A.scale(2.0);
    mult_dscal(2.0, A_copy);
    mult_dscal_to(2.0, A_copy2, B);
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < A.col(); j++) {
            EXPECT_DOUBLE_EQ(0.0, A(i, j));
            EXPECT_DOUBLE_EQ(0.0, A_copy(i, j));
            EXPECT_DOUBLE_EQ(0.0, B(i, j));
        }
    }

    // scale = positive
    A.randomize(0, 1);
    A_copy2 = A_copy = A;
    MatrixXd A_mxd = Matrix_to_MatrixXd(A);
    A_mxd *= 2.0;
    A.scale(2.0);
    mult_dscal(2.0, A_copy);
    mult_dscal_to(2.0, A_copy2, B);
    check_data_equality_with_EigenMatrix(A_mxd, A);
    check_data_equality_with_EigenMatrix(A_mxd, A_copy);
    check_data_equality_with_EigenMatrix(A_mxd, B);

    // scale = 0.0;
    A_copy2 = A_copy = A;
    A_mxd *= 0.0;
    A.scale(0.0);
    mult_dscal(0.0, A_copy);
    mult_dscal_to(0.0, A_copy2, B);
    check_data_equality_with_EigenMatrix(A_mxd, A);
    check_data_equality_with_EigenMatrix(A_mxd, A_copy);
    check_data_equality_with_EigenMatrix(A_mxd, B);

    A_copy2 = A_copy = A;
    A_mxd *= 1e-300;
    A.scale(1e-300);
    mult_dscal(1e-300, A_copy);
    mult_dscal_to(1e-300, A_copy2, B);
    check_data_equality_with_EigenMatrix(A_mxd, A);
    check_data_equality_with_EigenMatrix(A_mxd, A_copy);
    check_data_equality_with_EigenMatrix(A_mxd, B);
}
