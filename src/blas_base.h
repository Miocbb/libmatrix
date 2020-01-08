/**
 * file: blas_base.h
 *
 * declearation for blas library.
 */

#ifndef _MATRIX_BLAS_BASE_H_
#define _MATRIX_BLAS_BASE_H_

namespace matrix {

namespace blas {

static int izero[] = {0};
static int ione[] = {1};
static double dzero[] = {0.0};
static double done[] = {1.0};

extern "C" void dgemm_(const char *transa, const char *transb,
                       const int *m, const int *n, const int *k,
                       const double *alpha,
                       const double *a, const int *lda,
                       const double *b, const int *ldb,
                       const double *beta, double *c, const int *ldc);

extern "C" void dsyev_(const char * jobz, const char *uplo,
                       const int *N, double *A, const int *lda,
                       double *eig, double *work, const int *lwork, int *info);

extern "C" void dscal_(const int *N, const double *alpha, double *a,
                       const int *lda);
}

}

#endif
