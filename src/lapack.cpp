
#include <string>
#include <sstream>
#include "matrix.h"
#include "lapack.h"
#include "lapack_base.h"
#include "exception.h"

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
        throw exception::DimensionError("Cannot make a non-square matrix to be random orthogonal.");
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
        std::stringstream msg;
        msg << "QR factorization failed. " << "The " << -info
            << "-th argument had an illegal value\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    }
    // retrieve Q matrix from dgeqp3
    lapack::dorgqr_(&n, &n, &n, Q.data(), &n, tau, work, &lwork, &info);
    if (info < 0) {
        std::stringstream msg;
        msg << "Retrieve Q matrix failed. " << "The " << -info
            << "-th argument had an illegal value\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
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
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError("Cannot diagonalize a matrix that is not square.");
    } else if (A.row() > eig.size()) {
        string msg{"Fail to diagonalize a symmetric matrix: eigenvector size is too small."};
        throw exception::DimensionError(A.row(), eig.size(), msg);
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw matrix::exception::MatrixException("Unkown label to access a symmetric matrix data: label=" + uplo);
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
        throw matrix::exception::MatrixOperationError(__FUNCTION__, "convergence failure");
    } else if (info < 0) {
        std::stringstream msg;
        msg << "Fail to diagonalize a symmetric matrix: "
            << "the " << -info << "-th argument had an illegal value.\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
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
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError("Cannot invert a matrix that is not square.");
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
        std::stringstream msg;
        msg << "The " << -info << "-th arguments had an illegal value.\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    } else if (info > 0) {
        std::stringstream msg;
        msg << "U(" << -info << "," << -info << ") is exactly zero; the matrix is"
            << " singular and its inverse could not be computed.\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
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
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError("Cannot invert a matrix that is not squared.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw matrix::exception::MatrixException("Unkown label to access a symmetric matrix data: label=" + uplo);
    }

    int n = A.row();
    int info = 0;
    // call LAPACK to invert the matrix
    lapack::dpotrf_(used_uplo.c_str(), &n, A.data(), &n, &info);
    lapack::dpotri_(used_uplo.c_str(), &n, A.data(), &n, &info);
    if (info < 0) {
        std::stringstream msg;
        msg << "The " << -info << "-th arguments had an illegal value.";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    } else if (info > 0) {
        std::stringstream msg;
        msg << "The (" << -info << "," << -info << ") element of the factor U is zero,"
            << " and the inverse cannot be computed.";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
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
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError("Cannot invert a matrix that is not squared.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw exception::MatrixException("Unknown label to access a symmetric matrix data: label=" + uplo);
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
        std::stringstream msg;
        msg << "The " << -info << "-th arguments had an illegal value.";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    } else if (info > 0) {
        std::stringstream msg;
        msg << "D(" << info << "," << info << ") = zero; the matrix is"
            << "singular and its inverse could not be computed.";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
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
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError("Cannot invert a matrix that is not squared.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw exception::MatrixException("Unknown label to access a symmetric matrix data: label=" + uplo);
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
        std::stringstream msg;
        msg << "The " << -info << "-th arguments had an illegal value.";
        throw exception::MatrixOperationError(__FUNCTION__, msg.str());
    } else if (info > 0) {
        std::stringstream msg;
        msg << "D(" << info << "," << info << ") = zero; the matrix is"
            << "singular and its inverse could not be computed.";
        throw exception::MatrixOperationError(__FUNCTION__, msg.str());
    }
    // make inverse matrix full.
    A.to_symmetric(uplo);
    return 0;
}

}
