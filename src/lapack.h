#ifndef _MATRIX_SRC_LAPACK_H_
#define _MATRIX_SRC_LAPACK_H_

#include "matrix.h"

namespace matrix {

/**
 * set the input matrix be an random orthogonal matrix by using QR factorization
 * with column pivoting.
 *
 * The lapack subroutine dgeqp3 is used.
 *
 * @ param[in, out] Q the input matrix. On exit, it stores a random orthogonal matrix.
 * @ param[in] using_fixed_seed if or not using the fixed seed to feed the random number
 *  generator. If it is true, the function behavior can be repeatble at different running time,
 *  otherwise, it will not.
 * @ return int 0 refers to success.
 */
int set_matrix_random_orthogonal(Matrix & Q, bool using_fixed_seed = true);

/**
 * wrapper of lapack dsyev function to diagonalize a symmetric matrix.
 *
 * the eigenvalues and eigenvectors of the input matrix are computed,
 * and the input matrix data is destroyed on exit.
 *
 * @ param [in] uplo "U" upper triangular and "L" for lower triangular.
 * @ param [in, out] A the matrix to be diagonalized. On exit, the matrix stores
 *  the orthonormal eigenvalues of the matrix A when succeed.
 * @ param [out] the eigenvalues in ascending order when succeed.
 */
int diagonalize_sym_matrix_dsyev(const string & uplo, Matrix & A, vector<double> & eig);

/**
 * Invert a general matrix based on LAPACK dgetri, which is based on
 * LU factorization computed by LAPACK dgetrf.
 *
 * @ param[in,out] A the input general matrix. On exit, if succeed, it stores the inverse of
 *  the original matrix A.
 * @ return integer 0 refers to success.
 */
int invert_gen_matrix_dgetri(Matrix & A);

/**
 * Invert a real symmetric positive definite (spd) matrix by LAPACK dpotri,
 * which is based on Cholesky factorization computed by LAPACK dpotrf.
 *
 * The positive definite propoty of the input matrix is not checked when calling this function!
 *
 * @ param[in] uplo "U" using the upper triangular part; "L" using the lower triangular part.
 * @ param[in,out] A the input spd matrix. On exit, if succeed, it stores the inverse of
 *  the original matrix A.
 * @ return integer 0 refers to success.
 */
int invert_spd_matrix_dpotri(const string & uplo, Matrix & A);

/**
 * Invert a real symmetric indefinite matrix by LAPACK dsytri,
 * which is based on factorization A = U*D*U**T or A = L*D*L**T
 * computed by LAPACK dsytrf.
 *
 * @ param[in] uplo "U" using the upper triangular part; "L" using the lower triangular part.
 * @ param[in,out] A the input spd matrix. On exit, if succeed, it stores the inverse of
 *  the original matrix A.
 * @ return integer 0 refers to success.
 */
int invert_sym_matrix_dsytri(const string & uplo, Matrix &A);

/**
 * Invert a real symmetric indefinite matrix by LAPACK dsytri_rook,
 * which is based on factorization A = U*D*U**T or A = L*D*L**T
 * computed by LAPACK dsytrf_rook.
 *
 * @ param[in] uplo "U" using the upper triangular part; "L" using the lower triangular part.
 * @ param[in,out] A the input spd matrix. On exit, if succeed, it stores the inverse of
 *  the original matrix A.
 * @ return integer 0 refers to success.
 */
int invert_sym_matrix_dsytri_rook(const string & uplo, Matrix &A);

}   // namespace matrix

#endif  // _MATRIX_SRC_LAPACK_H_H
