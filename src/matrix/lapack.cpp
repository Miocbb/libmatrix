#include <matrix/details/exception.h>
#include <matrix/details/lapack.h>
#include <matrix/details/matrix.h>
#include <memory>
#include <sstream>
#include <string>

#include "lapack_base.h"

namespace matrix {

/**
 * @details Set the input matrix be an random orthogonal matrix by using QR
 * factorization with column pivoting. The lapack subroutine `dgeqp3` is used.
 */
int set_matrix_random_orthogonal(Matrix &Q, bool using_fixed_seed)
{
    if (!Q.is_square()) {
        throw exception::DimensionError(
            "Cannot make a non-square matrix to be random orthogonal.");
    }

    // first generate a random matrix.
    if (using_fixed_seed)
        Q.randomize_seed_fixed(0, 1);
    else
        Q.randomize(0, 1);

    // QR factorization to get Q matrix.
    int n = Q.col();
    int lwork = -1;
    int info = 0;
    double work_opt = 0;
    vector<int> jpvt(n);
    vector<double> tau(n);
    // query work space.
    lapack::dgeqp3_(&n, &n, Q.data(), &n, jpvt.data(), tau.data(), &work_opt,
                    &lwork, &info);
    lwork = (int)work_opt;
    vector<double> work(lwork);
    // do qr factorization.
    lapack::dgeqp3_(&n, &n, Q.data(), &n, jpvt.data(), tau.data(), work.data(),
                    &lwork, &info);
    if (info < 0) {
        std::stringstream msg;
        msg << "QR factorization failed. "
            << "The " << -info << "-th argument had an illegal value\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    }
    // retrieve Q matrix from dgeqp3
    lapack::dorgqr_(&n, &n, &n, Q.data(), &n, tau.data(), work.data(), &lwork,
                    &info);
    if (info < 0) {
        std::stringstream msg;
        msg << "Retrieve Q matrix failed. "
            << "The " << -info << "-th argument had an illegal value\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    }

    return 0;
}

int diagonalize_sym_matrix_dsyev(const string &uplo, Matrix &A,
                                 vector<double> &eig)
{
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError(
            "Cannot diagonalize a matrix that is not square.");
    } else if (A.row() > eig.size()) {
        string msg{"Fail to diagonalize a symmetric matrix: eigenvector size "
                   "is too small."};
        throw exception::DimensionError(A.row(), eig.size(), msg);
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw matrix::exception::MatrixException(
            "Unkown label to access a symmetric matrix data: label=" + uplo);
    }

    int row = A.row();
    int n = row;
    int lda = row;
    int info = 0;
    double wkopt = 0.0;
    // Query and allocate the optimal workspace
    int lwork = -1;
    lapack::dsyev_("V", used_uplo.c_str(), &n, A.data(), &lda, eig.data(),
                   &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    vector<double> work(lwork);

    // Solve eigenvalue decomposition.
    lapack::dsyev_("V", used_uplo.c_str(), &n, A.data(), &lda, eig.data(),
                   work.data(), &lwork, &info);

    // Check exit status
    if (info > 0) {
        throw exception::MatrixOperationError(__FUNCTION__,
                                              "convergence failure");
    } else if (info < 0) {
        std::stringstream msg;
        msg << "Fail to diagonalize a symmetric matrix: "
            << "the " << -info << "-th argument had an illegal value.\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    }

    return 0;
}

int invert_gen_matrix_dgetri(Matrix &A)
{
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError(
            "Cannot invert a matrix that is not square.");
    }

    int n = A.row();
    int lwork = n;
    int info = 0;
    vector<double> work(lwork);
    vector<int> ipiv(n);
    lapack::dgetrf_(&n, &n, A.data(), &n, ipiv.data(), &info);
    lapack::dgetri_(&n, A.data(), &n, ipiv.data(), work.data(), &lwork, &info);

    if (info < 0) {
        std::stringstream msg;
        msg << "The " << -info << "-th arguments had an illegal value.\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    } else if (info > 0) {
        std::stringstream msg;
        msg << "U(" << -info << "," << -info
            << ") is exactly zero; the matrix is"
            << " singular and its inverse could not be computed.\n";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    }
}

/**
 * @note The positive definite propoty (spd) of the input matrix is not checked
 * when calling this function!
 */
int invert_spd_matrix_dpotri(const string &uplo, Matrix &A)
{
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError(
            "Cannot invert a matrix that is not squared.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw matrix::exception::MatrixException(
            "Unkown label to access a symmetric matrix data: label=" + uplo);
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
        msg << "The (" << -info << "," << -info
            << ") element of the factor U is zero,"
            << " and the inverse cannot be computed.";
        throw matrix::exception::MatrixOperationError(__FUNCTION__, msg.str());
    }
    // make inverse matrix full. `dpotri` only update half of the matrix.
    A.to_symmetric(uplo);
    return 0;
}

int invert_sym_matrix_dsytri(const string &uplo, Matrix &A)
{
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError(
            "Cannot invert a matrix that is not squared.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw exception::MatrixException(
            "Unknown label to access a symmetric matrix data: label=" + uplo);
    }

    int n = A.row();
    int info = 0;
    double workopt = 0;
    vector<int> ipiv(n);
    // query and allocate the optimal workspace.
    int lwork = -1;
    lapack::dsytrf_(used_uplo.c_str(), &n, A.data(), &n, ipiv.data(), &workopt,
                    &lwork, &info);
    lwork = (int)workopt;
    vector<double> work(lwork);
    // call LAPACK to invert the matrix
    lapack::dsytrf_(used_uplo.c_str(), &n, A.data(), &n, ipiv.data(),
                    work.data(), &lwork, &info);
    lapack::dsytri_(used_uplo.c_str(), &n, A.data(), &n, ipiv.data(),
                    work.data(), &info);
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

int invert_sym_matrix_dsytri_rook(const string &uplo, Matrix &A)
{
    if (A.size() == 0) {
        return 0;
    } else if (!A.is_square()) {
        throw exception::DimensionError(
            "Cannot invert a matrix that is not squared.");
    }

    string used_uplo;
    if (uplo == "U") {
        used_uplo = "L";
    } else if (uplo == "L") {
        used_uplo = "U";
    } else {
        throw exception::MatrixException(
            "Unknown label to access a symmetric matrix data: label=" + uplo);
    }

    int n = A.row();
    int info = 0;
    double workopt = 0;
    vector<int> ipiv(n);
    // query and allocate the optimal workspace.
    int lwork = -1;
    lapack::dsytrf_rook_(used_uplo.c_str(), &n, A.data(), &n, ipiv.data(),
                         &workopt, &lwork, &info);
    lwork = (int)workopt;
    vector<double> work(lwork);
    // call LAPACK to invert the matrix
    lapack::dsytrf_rook_(used_uplo.c_str(), &n, A.data(), &n, ipiv.data(),
                         work.data(), &lwork, &info);
    lapack::dsytri_rook_(used_uplo.c_str(), &n, A.data(), &n, ipiv.data(),
                         work.data(), &info);
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

} // namespace matrix
