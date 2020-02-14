#ifndef _MATRIX_BLAS_H_
#define _MATRIX_BLAS_H_

#include "matrix.h"

namespace matrix {

using std::string;

/**
 * wrapper of blas dgemm function for general matrix multiplication.
 *
 * calculate C = alpha * op(A) * op(B) + beta * C.
 * op(A) represents the matrix with an operation acted on it.
 * op(A) = A or op(A) = A^T.
 *
 * @ param[in] alpha scalar coefficient on op(A) * op(B).
 * @ param[in] A matrix view that represents op(A).
 * @ param[in] op_A operation acting on matrix A.
 * @ param[in] B matrix view that represents op(B).
 * @ param[in] op_B operation acting on matrix B.
 * @ param[in] beta scalar coefficient on matrix C.
 * @ param[out] C matrix C.
 */
int mult_dgemm(const double alpha, const Matrix & A, const string & op_A,
               const Matrix & B, const string & op_B, const double beta,
               Matrix & C);

/**
 * calculate C = A * B * A^T.
 */
int mult_dgemm_ABAT(const Matrix &A, const Matrix &B, Matrix &C);

/**
 * wrapper of blas dscal function to scale matrix by a constant.
 *
 * @ param[in] alpha the scalar coefficient.
 * @ param[in] A the matrix to be scaled.
 * @ return int 0 for success others for failure.
 */
inline int mult_dscal(const double alpha, Matrix & A)
{
    A.scale(alpha);
    return 0;
}

/**
 * scale a matrix to another matrix by a constant.
 *
 * @ param[in] alpha the scalar coefficient.
 * @ param[in] A the matrix use for scaling.
 * @ param[out] B the scaled matrix.
 * @ return int 0 for success others for failure.
 */
int mult_dscal_to(const double alpha, const Matrix & A, Matrix & B);

}

#endif
