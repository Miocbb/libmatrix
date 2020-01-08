
#ifndef _MATRIX_LAPACK_H_
#define _MATRIX_LAPACK_H_

#include "matrix.h"

namespace matrix {

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
int diag_sym_matrix(const string & uplo, Matrix & A, vector<double> & eig);

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

}

#endif
