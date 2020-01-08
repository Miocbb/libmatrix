/**
 * file: lapack_base.h
 *
 * declearation for lapack library.
 */

#ifndef _MATRIX_LAPACK_BASE_H_
#define _MATRIX_LAPACK_BASE_H_

namespace matrix {

namespace lapack {

extern "C" void dsyev_(const char * jobz, const char *uplo,
                       const int *N, double *A, const int *lda,
                       double *eig, double *work, const int *lwork, int *info);

}

}

#endif
