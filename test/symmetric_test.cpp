#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

/**
 * Test matrix symmetrization method.
 */
TEST(ToSymmetricTest, test)
{
    size_t row = 3;
    size_t col = 3;
    Matrix A(row, col);
    ASSERT_EQ(row, col);

    // use lower part to symmetrize.
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < i; j++) {
            A(i, j) = i;
        }
    }
    Matrix B = A;
    A.symmetrize_lower_to_upper();
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < i; j++) {
            EXPECT_DOUBLE_EQ(A(i, j), A(j, i)) << "Use lower part: not symmetric at position "
                << i << "," << j << "]\n";
            EXPECT_DOUBLE_EQ(A(i, j), B(i, j)) << "Use lower part: wrong value at position "
                << i << "," << j << "]\n";
        }
    }

    // use upper part to symmetrize.
    for (size_t i = 0; i < row; i++) {
        for (size_t j = i; j < col; j++) {
            A(i, j) = i * 10;
        }
    }
    B = A;
    A.symmetrize_upper_to_lower();
    for (size_t i = 0; i < row; i++) {
        for (size_t j = i; j < col; j++) {
            EXPECT_DOUBLE_EQ(A(i, j), A(j, i)) << "Use upper part: not symmetric at position "
                << i << "," << j << "]\n";
            EXPECT_DOUBLE_EQ(A(i, j), B(i, j)) << "Use upper part: wrong value at position "
                << i << "," << j << "]\n";
        }
    }
}
