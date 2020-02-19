#include "utils.h"
#include <algorithm>
#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <random>
#include <vector>

using std::vector;
using matrix::Matrix;

void verify(Matrix &A, Matrix &Q_ref, vector<double> &eig_ref, Matrix &Q_calc,
            vector<double> &eig_calc)
{
    const size_t n = A.row();
    // verify that eigenvector matrix. Check (Q_calc * Q_ref)^2 == I
    Matrix eigM_Q(n, n);
    matrix::mult_dgemm(1.0, Q_ref, "N", Q_calc, "T", 0.0, eigM_Q);
    Matrix eigM_Q_2(3, 3);
    matrix::mult_dgemm(1.0, eigM_Q, "N", eigM_Q, "T", 0.0, eigM_Q_2);
    // EXPECT_TRUE(eigM_Q.is_identity());
    auto status1 = eigM_Q_2.is_identity();
    EXPECT_TRUE(status1);
    if (!status1) {
        std::cout << "Q_ref:\n";
        Q_ref.show_full();
        std::cout << "Q_calc:\n";
        Q_calc.show_full();
    }

    // verify eigenvalues.
    Matrix eigV_calc(1, n, eig_calc);
    Matrix eigV_ref(1, n, eig_ref);
    auto status2 = eigV_calc.is_equal_to(eigV_ref);
    EXPECT_TRUE(status2);
    if (!status2) {
        std::cout << "eigenvalue ref:\n";
        eigV_ref.show_full();
        std::cout << "eigenvalue calc:\n";
        eigV_calc.show_full();
    }

    // verify calculated EVD can restore the original input matrix.
    Matrix A2(n, n);
    Matrix D(n, n);
    for (size_t i = 0; i < D.size(); ++i)
        D(i, i) = eig_calc[i];
    matrix::mult_dgemm_ATBA(Q_calc, D, A2);
    EXPECT_TRUE(A.is_equal_to(A2));
}

TEST(DiagonalizeTest, dsyev_test)
{
    Matrix Q(3, 3);
    matrix::set_matrix_random_orthogonal(Q);
    Matrix DiagM(3, 3);

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<> dist6(-10, 10);
    std::vector<double> rand_vec;
    for (size_t i = 0; i < Q.row(); ++i)
        rand_vec.push_back(dist6(rng));
    std::sort(rand_vec.begin(), rand_vec.end());
    for (size_t i = 0; i < Q.row(); ++i)
        DiagM(i, i) = rand_vec[i];

    // make symmetric matrix that can be diagonalized.
    Matrix A(3, 3);
    std::vector<double> eigV(3);
    matrix::mult_dgemm_ATBA(Q, DiagM, A);

    // do diagonalization with Upper part.
    Matrix Q_calc = A;
    matrix::diagonalize_sym_matrix_dsyev("U", Q_calc, eigV);
    verify(A, Q, rand_vec, Q_calc, eigV);

    // do diagonalization with Lower part.
    Q_calc = A;
    matrix::diagonalize_sym_matrix_dsyev("L", Q_calc, eigV);
    verify(A, Q, rand_vec, Q_calc, eigV);
}
