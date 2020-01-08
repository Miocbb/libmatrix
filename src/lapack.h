
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
int eig_dsyev(const string & uplo, Matrix & A, vector<double> & eig);

}

#endif
