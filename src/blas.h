#ifndef _MATRIX_BLAS_H_
#define _MATRIX_BLAS_H_

#include "matrix.h"
#include "io.h"

namespace matrix {

using std::string;
using matrix::Matrix;

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
}

#endif
