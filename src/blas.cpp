#include <string>
#include "matrix.h"
#include "blas.h"
#include "blas_base.h"
#include "io.h"

namespace matrix {

using std::string;

/**
 * General matrix multiplication: C = alpha * A * B + beta * C.
 *
 * @ param[in] alpha scalar coefficient on A * B.
 * @ param[in] A matrix A.
 * @ param[in] B matrix B.
 * @ param[in] beta scalar coefficient on matrix C.
 * @ param[out] C matrix C.
 */
static int mult_dgemm_NN(double alpha, const Matrix & A, const Matrix & B,
                         double beta, Matrix & C)
{
    // check dimension
    // A: M x N
    // B: N x K
    // AB: M x K
    // BT AT: K x M
    int M = A.row();
    int N = A.col();
    int K = B.col();
    // dimension check.
    if (N != B.row()) {
        // dimension check for A and B:
        sig_err("Error in matrix::mult_dgemm_NN(): dimension error between matrix A and B.\n");
    } else if (M != C.row() || K != C.col()) {
        // dimension check for AB and C:
        sig_err("Error in matrix::mult_dgemm_NN(): dimension error between matrix (AB) and C.\n");
    }
    // calculate (AB)^T=(B^T A^T) by dgemm to get (AB) stored in row-wise matrix.
    blas::dgemm_("N", "N", &K, &M, &N,
                 &alpha,
                 B.data(), &K,
                 A.data(), &N,
                 &beta, C.data(), &K);
    return 0;
}

/**
 * General matrix multiplication: C = alpha * A^T * B^T + beta * C.
 *
 * @ param[in] alpha scalar coefficient on A * B.
 * @ param[in] A matrix A.
 * @ param[in] B matrix B.
 * @ param[in] beta scalar coefficient on matrix C.
 * @ param[out] C matrix C.
 */
static int mult_dgemm_TT(double alpha, const Matrix & A, const Matrix & B,
                         double beta, Matrix & C)
{
    // check dimension
    // A: M x N
    // B: K x M
    // A^T B^T: N x K
    // BA: K x N
    int M = A.row();
    int N = A.col();
    int K = B.row();
    if (M != B.col()) {
        // dimension check for A and B:
        sig_err("Error in matrix::mult_dgemm_TT(): dimension error between matrix A^T and B^T.\n");
    } else if (N != C.row() || K != C.col()) {
        // dimension check for AB and C:
        sig_err("Error in matrix::mult_dgemm_TT(): dimension error between matrix (A^T B^T) and C.\n");
    }
    // calculate BA by dgemm to get (A^T B^T) stored in row-wise matrix.
    blas::dgemm_("T", "T", &K, &N, &M,
                 &alpha,
                 B.data(), &M,
                 A.data(), &N,
                 &beta, C.data(), &K);
    return 0;
}

/**
 * General matrix multiplication: C = alpha * A * B^T + beta * C.
 *
 * @ param[in] alpha scalar coefficient on A * B.
 * @ param[in] A matrix A.
 * @ param[in] B matrix B.
 * @ param[in] beta scalar coefficient on matrix C.
 * @ param[out] C matrix C.
 */
static int mult_dgemm_NT(double alpha, const Matrix & A, const Matrix & B,
                         double beta, Matrix & C)
{
    // check dimension
    // A: M x N
    // B: K x N
    // A B^T: M x K
    // B A^T: K x M
    int M = A.row();
    int N = A.col();
    int K = B.row();
    if (N != B.col()) {
        // dimension check for A and B:
        sig_err("Error in matrix::mult_dgemm_NT(): dimension error between matrix A and B^T.\n");
    } else if (M != C.row() || K != C.col()) {
        // dimension check for AB and C:
        sig_err("Error in matrix::mult_dgemm_NT(): dimension error between matrix (A B^T) and C.\n");
    }
    // calculate B A^T by dgemm to get (A B^T) stored in row-wise matrix.
    blas::dgemm_("T", "N", &K, &M, &N,
                 &alpha,
                 B.data(), &N,
                 A.data(), &N,
                 &beta, C.data(), &K);
    return 0;
}

/**
 * General matrix multiplication: C = alpha * A^T * B + beta * C.
 *
 * @ param[in] alpha scalar coefficient on A * B.
 * @ param[in] A matrix A.
 * @ param[in] B matrix B.
 * @ param[in] beta scalar coefficient on matrix C.
 * @ param[out] C matrix C.
 */
static int mult_dgemm_TN(double alpha, const Matrix & A, const Matrix & B,
                         double beta, Matrix & C)
{
    // A: M x N
    // B: M x K
    // A^T B: N x K
    // B^T A: K x N
    int M = A.row();
    int N = A.col();
    int K = B.col();
    if (M != B.row()) {
        // dimension check for A and B:
        sig_err("Error in matrix::mult_dgemm_TN(): dimension error between matrix A^T and B.\n");
    } else if (N != C.row() || K != C.col()) {
        // dimension check for AB and C:
        sig_err("Error in matrix::mult_dgemm_TN(): dimension error between matrix (A^T B) and C.\n");
    }
    // calculate B^T A by dgemm to get (A B^T) stored in row-wise matrix.
    blas::dgemm_("N", "T", &K, &N, &M,
                 &alpha,
                 B.data(), &K,
                 A.data(), &N,
                 &beta, C.data(), &K);
    return 0;
}

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
               Matrix & C)
{
    if (op_A == "N" && op_B == "N") {
        return mult_dgemm_NN(alpha, A, B, beta, C);
    } else if (op_A == "N" && op_B == "T") {
        return mult_dgemm_NT(alpha, A, B, beta, C);
    } else if (op_A == "T" && op_B == "N") {
        return mult_dgemm_TN(alpha, A, B, beta, C);
    } else if (op_A == "T" && op_B == "T") {
        return mult_dgemm_TT(alpha, A, B, beta, C);
    } else {
        sig_err("Error in matrix::mult_dgemm(): unknow operation on matrix\n");
    }
}

/**
 * scale a matrix to another matrix by a constant.
 *
 * @ param[in] alpha the scalar coefficient.
 * @ param[in] A the matrix use for scaling.
 * @ param[out] B the scaled matrix.
 * @ return int 0 for success others for failure.
 */
int mult_dscal_to(const double alpha, const Matrix & A, Matrix & B)
{
    if (A.row() != B.row() || A.col() != B.col()) {
        sig_err("Error in matrix::mult_dscal_to(), matrix dimension mismatched.\n");
    }
    if (alpha == 1.0) {
        B = A;
        return 0;
    } else if (alpha == 0.0) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < A.size(); i++) {
            B.data()[i] = 0.0;
        }
    } else {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < A.size(); i++) {
            B.data()[i] = A.data()[i] * alpha;
        }
    }
    return 0;
}

}
