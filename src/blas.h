/**
 * @file
 * @brief declaration of matrix library functions relate to functions in blas
 * library.
 */

#ifndef _MATRIX_SRC_BLAS_H_
#define _MATRIX_SRC_BLAS_H_

#include "matrix.h"

namespace matrix {

using std::string;

/**
 * @brief wrapper of blas `dgemm` function for general matrix multiplication.
 *
 * @par Purpose
 * calculate C = alpha * op(A) * op(B) + beta * C.\n
 * op(A) represents the matrix with an operation acted on it.\n
 * op(A) = A or op(A) = A^T.
 *
 * @param [in] alpha: scalar coefficient on op(A) * op(B).
 * @param [in] A: matrix view that represents op(A).
 * @param [in] op_A:  operation acting on matrix A.
 * @param [in] B: matrix view that represents op(B).
 * @param [in] op_B: operation acting on matrix B.
 * @param [in] beta: scalar coefficient on matrix C.
 * @param [out] C: matrix C.
 */
int mult_dgemm(const double alpha, const Matrix &A, const string &op_A,
               const Matrix &B, const string &op_B, const double beta,
               Matrix &C);

/**
 * @brief convenient function wrapper for three general matrix multiplication.
 *
 * @par Purpose
 * calculate C = A * B * A^T.
 *
 * @param [in] A: matrix A.
 * @param [in] B: matrix B.
 * @param [in] C: matrix C.
 */
int mult_dgemm_ABAT(const Matrix &A, const Matrix &B, Matrix &C);

/**
 * @brief convenient function wrapper for three general matrix multiplication.
 *
 * @par Purpose
 * calculate C = A^T * B * A.
 *
 * @param [in] A: matrix A.
 * @param [in] B: matrix B.
 * @param [in] C: matrix C.
 */
int mult_dgemm_ATBA(const Matrix &A, const Matrix &B, Matrix &C);

/**
 * @brief wrapper of blas dscal function to scale matrix by a constant.
 *
 * @par Purpose
 * calculate A = alpha * A
 *
 * @param [in] alpha the scalar coefficient.
 * @param [in, out] A the matrix to be scaled. On exit, matrix A is updated.
 * @return int 0 for success others for failure.
 */
inline int mult_dscal(const double alpha, Matrix &A)
{
    A.scale(alpha);
    return 0;
}

/**
 * @brief scale a matrix to another matrix by a constant.
 *
 * @par Purpose
 * calculate B = alpha * A
 *
 * @param [in] alpha the scalar coefficient.
 * @param [in] A the matrix use for scaling.
 * @param [in, out] B On exit, matrix B is overwritten and stores the scaled
 * matrix.
 * @return int 0 for success others for failure.
 */
int mult_dscal_to(const double alpha, const Matrix &A, Matrix &B);

} // namespace matrix

#endif // _MATRIX_SRC_BLAS_H_
