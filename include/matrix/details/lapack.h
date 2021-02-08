/**
 * @file lapack.h
 * @brief declaration of matrix library functions relate to lapack function.
 */
#ifndef _MATRIX_INCLUDE_MATRIX_DETAILS_LAPACK_H_
#define _MATRIX_INCLUDE_MATRIX_DETAILS_LAPACK_H_

#include "matrix.h"

namespace matrix {

/**
 * @brief Set the input matrix be an random orthogonal matrix.
 * @param [in, out] A: The input matrix. On exit, it stores a random orthogonal
 * matrix.
 * @param [in] using_fixed_seed: If or not using the fixed seed to feed the
 * random number generator. If it is true, the function behavior can be
 * repeatble at different running time, otherwise, it will not.
 * @return int: 0 for success, and others for failure.
 */
int set_matrix_random_orthogonal(Matrix &A, bool using_fixed_seed = true);

/**
 * @brief Wrapper of lapack `dsyev` function to diagonalize a symmetric matrix.
 *
 * @par Purpose
 * Calculate A = Q^T * D * Q, where D is a diagonal matrix, and Q is an
 * orthonormal matrix.
 *
 * @param [in] uplo: "U": only the upper triangular will be refereed.\n
 * "L": only the lower triangular will be refereed.
 * @param [in, out] A: The matrix to be diagonalized. On exit, the matrix stores
 * the orthonormal eigenvalues matrix Q.
 * @param [out] eig: The eigenvalues in ascending order when succeed.
 * @return int: 0 for success, and others for failure.
 *
 * @note On successful exit, matrix `A` stores the eigenvalues matrix Q,
 * that is each eigenvector stores continuously in memory.
 */
int diagonalize_sym_matrix_dsyev(const string &uplo, Matrix &A,
                                 vector<double> &eig);

/**
 * @brief Invert a general matrix based on lapack `dgetri`, which is based on
 * LU factorization computed by lapack `dgetrf`.
 *
 * @param[in,out] A: The input general matrix. On exit, if succeed, it stores
 * the inverse of the original matrix A.
 * @return int: 0 for success, and others for failure.
 */
int invert_gen_matrix_dgetri(Matrix &A);

/**
 * @brief Invert a real symmetric positive definite (spd) matrix by lapack
 * `dpotri`, which is based on Cholesky factorization computed by lapack
 * `dpotrf`.
 *
 * @param [in] uplo: "U": only the upper triangular will be refereed.\n
 *  "L": only the lower triangular will be refereed.
 * @param[in,out] A: The input spd matrix. On exit, if succeed, it stores the
 * inverse of the original matrix A.
 * @return int: 0 for success, and others for failure.
 */
int invert_spd_matrix_dpotri(const string &uplo, Matrix &A);

/**
 * @brief Invert a real symmetric indefinite matrix by lapack `dsytri`,
 * which is based on factorization A = U*D*U**T or A = L*D*L**T
 * computed by lapack `dsytrf`.
 *
 * @param [in] uplo: "U": only the upper triangular will be refereed.\n
 * "L": only the lower triangular will be refereed.
 * @param[in,out] A: The input spd matrix. On exit, if succeed, it stores the
 * inverse of the original matrix A.
 * @return int: 0 for success, and others for failure.
 */
int invert_sym_matrix_dsytri(const string &uplo, Matrix &A);

/**
 * @brief Invert a real symmetric indefinite matrix by lapack `dsytri_rook`,
 * which is based on factorization A = U*D*U**T or A = L*D*L**T
 * computed by lapack `dsytrf_rook`.
 *
 * @param [in] uplo: "U": only the upper triangular will be refereed.\n
 * "L": only the lower triangular will be refereed.
 * @param [in,out] A: The input matrix. On exit, if succeed, it stores the
 * inverse of the original matrix A.
 * @return int: 0 for success, and others for failure.
 */
int invert_sym_matrix_dsytri_rook(const string &uplo, Matrix &A);

} // namespace matrix

#endif // _MATRIX_SRC_LAPACK_H_H
