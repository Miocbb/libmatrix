#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;
using matrix::diagonalize_sym_matrix_dsyev;

/**
 * test symmetric matrix diagonalization function.
 */
TEST(DiagonalizeSymmetricMatrixTest, dsyev_test)
{
    /**
     * trivial case: diagonalize a diagonal matrix.
     * A = [0, 0,
     *      0, 1]
     */
    Matrix A(2, 2);
    vector<double> A_diag= {0, 1};
    for (size_t i = 0; i < 2; i++) {
        A(i, i) = A_diag[i];
    }

    Matrix A_eigvec;
    vector<double> A_eigval(2);

    // test upper triangular.
    A_eigvec = A;
    diagonalize_sym_matrix_dsyev("U", A_eigvec, A_eigval);
    for (size_t i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(A_eigval[i], A_diag[i]);
    }
    //A_eigvec.show_full();
    EXPECT_TRUE(A_eigvec.is_identity());     // eigenvector matrix should be identity.
    // test lower triangular.
    A_eigvec = A;
    diagonalize_sym_matrix_dsyev("L", A_eigvec, A_eigval);
    for (size_t i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(A_eigval[i], A_diag[i]);
    }
    EXPECT_TRUE(A_eigvec.is_identity());     // eigenvector matrix should be identity.

    /**
     * A = [1, 1,
     *      1, 1]
     */
    A.fill_all(1.0);
    A_diag = {0.0, 2.0};
    // test upper triangular.
    A_eigvec = A;
    diagonalize_sym_matrix_dsyev("U", A_eigvec, A_eigval);
    for (size_t i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(A_eigval[i], A_diag[i]);
    }
    // test lower triangular.
    A_eigvec = A;
    diagonalize_sym_matrix_dsyev("L", A_eigvec, A_eigval);
    for (size_t i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(A_eigval[i], A_diag[i]);
    }

    /**
     * A = [1, xxx,
     *      1, 1]
     */
    A.fill_all(1.0);
    A(0, 1) = 999;
    A_diag = {0.0, 2.0};
    // test upper triangular.
    A_eigvec = A;
    diagonalize_sym_matrix_dsyev("L", A_eigvec, A_eigval);
    for (size_t i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(A_eigval[i], A_diag[i]);
    }

    /**
     * A = [1, 1,
     *      xxx, 1]
     */
    A.fill_all(1.0);
    A(1, 0) = 999;
    A_diag = {0.0, 2.0};
    // test upper triangular.
    A_eigvec = A;
    diagonalize_sym_matrix_dsyev("U", A_eigvec, A_eigval);
    for (size_t i = 0; i < 2; i++) {
        EXPECT_DOUBLE_EQ(A_eigval[i], A_diag[i]);
    }
}
