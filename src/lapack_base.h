/**
 * file: lapack_base.h
 *
 * declearation for lapack library.
 */

#ifndef _MATRIX_SRC_LAPACK_BASE_H_
#define _MATRIX_SRC_LAPACK_BASE_H_

namespace matrix {

namespace lapack {

extern "C" void dsyev_(const char *jobz, const char *uplo,
                       const int *N, double *A, const int *lda,
                       double *eig, double *work, const int *lwork, int *info);
extern "C" void dgeqp3_(const int *M, const int *N, double *A, const int *lda,
                        int *jpvt, double *tau, double *work, int *lwork, int *info);
extern "C" void dorgqr_(const int *M, const int *N, const int *K,
                        double *A, const int *lda, double *tau, double *work, int *lwork, int *info);
extern "C" void dgetrf_(const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info);
extern "C" void dgetri_(const int *n, double *a, const int *lda, int *ipiv, double *work, int *lwork, int *info);
extern "C" void dpotrf_(const char *uplo, const int *n, double *a, const int *lda, int *info);
extern "C" void dpotri_(const char *uplo, const int *n, double *a, const int *lda, int *info);
extern "C" void dsytrf_(const char *uplo, const int *n, double *a, const int *lda, int *ipiv, double *work, int *lwork, int*info);
extern "C" void dsytri_(const char *uplo, const int *n, double *a, const int *lda, int *ipiv, double *work, int*info);
extern "C" void dsytrf_rook_(const char *uplo, const int *n, double *a, const int *lda, int *ipiv, double *work, int *lwork, int*info);
extern "C" void dsytri_rook_(const char *uplo, const int *n, double *a, const int *lda, int *ipiv, double *work, int*info);

}   // namespace matrix::lapack
}   // namespace matrix

#endif  // _MATRIX_SRC_LAPACK_BASE_H_
