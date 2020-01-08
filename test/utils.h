#ifndef _MATRIX_TEST_UTILS_H_
#define _MATRIX_TEST_UTILS_H_
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <vector>
#include <matrix/matrix.h>

using matrix::Matrix;
using matrix::mult_dgemm;
using Eigen::MatrixXd;

inline MatrixXd Matrix_to_MatrixXd(Matrix & m)
{
    MatrixXd mxd(m.row(), m.col());
    for (size_t i = 0; i < m.row(); i++) {
        for (size_t j = 0; j < m.col(); j++) {
            mxd(i, j) = m(i, j);
        }
    }
    return mxd;
}

inline Matrix MatrixXd_to_Matrix(MatrixXd & mxd)
{
    Matrix mat(mxd.rows(), mxd.cols());
    for (size_t i = 0; i < mxd.rows(); i++) {
        for (size_t j = 0; j < mxd.cols(); j++) {
            mat(i, j) = mxd(i, j);
        }
    }
    return mat;
}

inline void check_data_equality_with_EigenMatrix(MatrixXd & ref, Matrix & rst)
{
    ASSERT_EQ(ref.rows(), rst.row()) << "can not compare matrix: row dimension mismatch.";
    ASSERT_EQ(ref.cols(), rst.col()) << "can not compare matrix: col dimension mismatch.";
    for (size_t i = 0; i < rst.row(); i++) {
        for (size_t j = 0; j < rst.col(); j++) {
            EXPECT_DOUBLE_EQ(ref(i, j), rst(i, j)) << "wrong element value at position ["
                << i << "," << j << "].\n";
        }
    }
}

#endif
