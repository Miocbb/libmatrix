
#include <string>
#include "matrix.h"
#include "lapack.h"
#include "lapack_base.h"

namespace matrix {

/**
 * set the input matrix be an random orthogonal matrix by using QR factorization
 * with column pivoting.
 *
 * The lapack subroutine dgeqp3 is used.
 *
 * @ param[in, out] Q the input matrix. On exit, it stores a random orthogonal matrix.
 * @ return int 0 refers to success.
 */
int set_matrix_random_orthogonal(Matrix & Q, bool using_fixed_seed)
{
    if (!Q.is_square()) {
        sig_err("Error in set_matrix_random_orthogonal: matrix is not a square matrix.\n");
    }

    // first generate a random matrix.
    if (using_fixed_seed) {
        Q.randomize_seed_fixed(0, 1);
    } else {
        Q.randomize(0, 1);
    }

    // QR factorization to get Q matrix.
    int n = Q.col();
    int *jpvt = new int[n];
    double *tau = new double[n];
    double work_opt = 0;
    int lwork = -1;
    double *work = nullptr;
    int info = 0;
    // query work space.
    lapack::dgeqp3_(&n, &n, Q.data(), &n, jpvt, tau, &work_opt, &lwork, &info);
    lwork = (int) work_opt;
    work = new double[lwork];
    // do qr factorization.
    lapack::dgeqp3_(&n, &n, Q.data(), &n, jpvt, tau, work, &lwork, &info);
    if (info < 0) {
        printf("Error in mtx_set_matrix_random_orthogonal:"
               "QR factorization failed."
               "The %d-th argument had an illegal value\n", -info);
        std::exit(EXIT_FAILURE);
    }
    // retrieve Q matrix from dgeqp3
    lapack::dorgqr_(&n, &n, &n, Q.data(), &n, tau, work, &lwork, &info);
    if (info < 0) {
        printf("Error in mtx_set_matrix_random_orthogonal:"
               "Retrieve Q matrix failed."
               "The %d-th argument had an illegal value\n", -info);
        std::exit(EXIT_FAILURE);
    }

    delete [] jpvt;
    delete [] tau;
    delete [] work;

    return 0;
}

/**
 * wrapper of blas dsyev function to diagonalize a symmetric matrix.
 *
 * the eigenvalues and eigenvectors of the input matrix are computed,
 * and the input matrix data is destroyed on exit.
 *
 * @ param [in] uplo "U" upper triangular and "L" for lower triangular.
 * @ param [in, out] A the matrix to be diagonalized. On exit, the matrix stores
 *  the orthonormal eigenvalues of the matrix A when succeed.
 * @ param [out] the eigenvalues in ascending order when succeed.
 */
int diagonalize_sym_matrix_dsyev(const string & uplo, Matrix & A, vector<double> & eig)
{
    if (!A.is_square()) {
        sig_err("Error to diagonalize a symmetric matrix: it is not even squared.");
    } else if (A.row() > eig.size()) {
        sig_err("Error to diagonalize a symmetric matrix: eigenvalue vector size is smaller than matrix dimension.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        sig_err("Error to diagonalize a symmetric matrix: unkown label to access matrix data.");
    }

    int row = A.row();
    int n = row;
    int lda = row;
    int info = 0;
    int lwork = -1;
    double wkopt = 0.0;
    double *w = eig.data();
    double *a = A.data();
    // Query and allocate the optimal workspace
    matrix::lapack::dsyev_("V", used_uplo.c_str(), &n, a, &lda, w, &wkopt, &lwork, &info);
    lwork = (int) wkopt;
    double *work = new double [lwork];

    // Solve eigenproblem
    matrix::lapack::dsyev_("V", used_uplo.c_str(), &n, a, &lda, w, work, &lwork, &info);
    delete [] work;

    // Check exit status
    if (info > 0) {
        sig_err("Error to diagonalize a symmetric matrix: diagonalization failed to converge.\n");
    } else if (info < 0) {
        std::cout << "Error to diagonalize a symmetric matrix:\n";
        std::cout << "the " << info << "-th argument had an illegal value.\n";
        std::exit(EXIT_FAILURE);
    }

    return 0;
}

/**
 * Invert a general matrix based on LAPACK dgetri, which is based on
 * LU factorization computed by LAPACK dgetrf.
 *
 * @ param[in,out] A the input general matrix. On exit, if succeed, it stores the inverse of
 *  the original matrix A.
 * @ return integer 0 refers to success.
 */
int invert_gen_matrix_dgetri(Matrix & A)
{
    if (!A.is_square()) {
        std::cout << "Error to invert a matrix: matrix is not square.";
        std::exit(EXIT_FAILURE);
    } else if (A.size() == 0) {
        std::cout << "Error to invert a matrix: matrix is empty.";
        std::exit(EXIT_FAILURE);
    }
    int n = A.row();
    int lwork = n;
    double *work = nullptr;
    work = new double[lwork];
    int *ipiv = nullptr;
    int info = 0;
    ipiv = new int[n];
    lapack::dgetrf_(&n, &n, A.data(), &n, ipiv, &info);
    lapack::dgetri_(&n, A.data(), &n, ipiv, work, &lwork, &info);
    delete [] work;
    delete [] ipiv;

    if (info < 0) {
        printf("Error in losc_ri_inverse_by_dgetri():\n"
               "the %d-th arguments had an illegal value.\n", -info);
        std::exit(EXIT_FAILURE);
    } else if (info > 0) {
        printf("Error in losc_ri_inverse_by_dgetri():\n"
               "U(%d,%d) is exactly zero; the matrix is"
               "singular and its inverse could not be computed.\n", info, info);
        std::exit(EXIT_FAILURE);
    }
}

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
int invert_spd_matrix_dpotri(const string & uplo, Matrix & A)
{
    if (!A.is_square()) {
        std::cout << "Error to invert a matrix: matrix is not square.";
        std::exit(EXIT_FAILURE);
    } else if (A.size() == 0) {
        std::cout << "Error to invert a matrix: matrix is empty.";
        std::exit(EXIT_FAILURE);
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        std::cout << "Error to invert a symmetric matrix: unkown label to access matrix data parts.";
        std::exit(EXIT_FAILURE);
    }

    int n = A.row();
    int info = 0;
    // call LAPACK to invert the matrix
    lapack::dpotrf_(used_uplo.c_str(), &n, A.data(), &n, &info);
    lapack::dpotri_(used_uplo.c_str(), &n, A.data(), &n, &info);
    if (info < 0) {
        printf("Error in invert_spd_matrix_dpotri():\n"
               "the %d-th arguments had an illegal value.\n", -info);
        std::exit(EXIT_FAILURE);
    } else if (info > 0) {
        printf("Error in invert_spd_matrix_dpotri():\n"
               "the (%d, %d) element of the factor U is zero, and the inverse could"
               "not be computed.", info, info);
        std::exit(EXIT_FAILURE);
    }
    // make inverse matrix full. `dpotri` only update half of the matrix.
    A.to_symmetric(uplo);
    return 0;
}

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
int invert_sym_matrix_dsytri(const string & uplo, Matrix &A)
{
    if (!A.is_square()) {
        std::cout << "Error to invert a matrix: matrix is not square.";
        std::exit(EXIT_FAILURE);
    } else if (A.size() == 0) {
        std::cout << "Error to invert a matrix: matrix is empty.";
        std::exit(EXIT_FAILURE);
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        std::cout << "Error to invert a symmetric matrix: unkown label to access matrix data parts.";
        std::exit(EXIT_FAILURE);
    }

    int n = A.row();
    int lwork = 0;
    double workopt = 0;
    double *work = nullptr;
    int *ipiv = new int[n];
    int info = 0;
    lwork = -1;
    // query and allocate the optimal workspace.
    lapack::dsytrf_(used_uplo.c_str(), &n, A.data(), &n, ipiv, &workopt, &lwork, &info);
    lwork = (int) workopt;
    work = new double[lwork];
    // call LAPACK to invert the matrix
    lapack::dsytrf_(used_uplo.c_str(), &n, A.data(), &n, ipiv, work, &lwork, &info);
    lapack::dsytri_(used_uplo.c_str(), &n, A.data(), &n, ipiv, work, &info);
    delete [] work;
    delete [] ipiv;
    if (info < 0) {
        printf("Error in invert_sym_matrix_dsytri():\n"
               "the %d-th arguments had an illegal value.\n", -info);
        std::exit(EXIT_FAILURE);
    } else if (info > 0) {
        printf("Error in invert_sym_matrix_dsytri():\n"
               "D(%d,%d) = zero; the matrix is"
               "singular and its inverse could not be computed.\n", info, info);
        std::exit(EXIT_FAILURE);
    }
    // make inverse matrix full.
    A.to_symmetric(uplo);
    return 0;
}

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
int invert_sym_matrix_dsytri_rook(const string & uplo, Matrix &A)
{
    if (!A.is_square()) {
        std::cout << "Error to invert a matrix: matrix is not square.";
        std::exit(EXIT_FAILURE);
    } else if (A.size() == 0) {
        std::cout << "Error to invert a matrix: matrix is empty.";
        std::exit(EXIT_FAILURE);
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        std::cout << "Error to invert a symmetric matrix: unkown label to access matrix data parts.";
        std::exit(EXIT_FAILURE);
    }

    int n = A.row();
    int lwork = 0;
    double workopt = 0;
    double *work = nullptr;
    int *ipiv = new int[n];
    int info = 0;
    lwork = -1;
    // query and allocate the optimal workspace.
    lapack::dsytrf_rook_(used_uplo.c_str(), &n, A.data(), &n, ipiv, &workopt, &lwork, &info);
    lwork = (int) workopt;
    work = new double[lwork];
    // call LAPACK to invert the matrix
    lapack::dsytrf_rook_(used_uplo.c_str(), &n, A.data(), &n, ipiv, work, &lwork, &info);
    lapack::dsytri_rook_(used_uplo.c_str(), &n, A.data(), &n, ipiv, work, &info);
    delete [] work;
    delete [] ipiv;
    if (info < 0) {
        printf("Error in invert_sym_matrix_dsytri_rook():\n"
               "the %d-th arguments had an illegal value.\n", -info);
        std::exit(EXIT_FAILURE);
    } else if (info > 0) {
        printf("Error in invert_sym_matrix_dsytri_rook():\n"
               "D(%d,%d) = zero; the matrix is"
               "singular and its inverse could not be computed.\n", info, info);
        std::exit(EXIT_FAILURE);
    }
    // make inverse matrix full.
    A.to_symmetric(uplo);
    return 0;
}

}
