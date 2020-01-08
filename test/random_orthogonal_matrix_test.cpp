#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;
using matrix::set_matrix_random_orthogonal;

/**
 * test random orthogonal matrix generation:
 * the function is matrix::set_matrix_random_orthogonal;
 */
TEST(SetRandomOrthogonalMatrixTest, check_orthogonality)
{
    Matrix Q(2, 2);
    set_matrix_random_orthogonal(Q, false);
    Matrix I(2, 2);
    matrix::mult_dgemm(1.0, Q, "N", Q, "T", 0.0, I);
    EXPECT_TRUE(I.is_identity());

    matrix::mult_dgemm(1.0, Q, "T", Q, "N", 0.0, I);
    EXPECT_TRUE(I.is_identity());

    // eye-ball the randomness.
    std::cout << "Below is a random orthogonal matrix. The random seed is not fixed. Does it look random?\n";
    Q.show_full();
}
