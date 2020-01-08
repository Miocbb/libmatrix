
#include <string>
#include "matrix.h"
#include "lapack.h"
#include "lapack_base.h"

namespace matrix {

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
int eig_dsyev(const string & uplo, Matrix & A, vector<double> & eig)
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

}
