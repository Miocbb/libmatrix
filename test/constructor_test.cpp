#include <gtest/gtest.h>
#include <matrix/matrix.h>

using matrix::Matrix;

/**
 * Test the only one constructor, which should construct a zero matrix.
 */
TEST(MatrixConstructorTest, initialization_test)
{
    // explicit constructor.
    size_t row = 2;
    size_t col = 3;
    Matrix A(row, col);
    EXPECT_EQ(row, A.row());
    EXPECT_EQ(col, A.col());
    EXPECT_EQ(row * col, A.size());
    // check zero initialization.
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < A.col(); j++) {
            EXPECT_DOUBLE_EQ(0.0, A(i, j)) << "Matrix zero initialization not working "
                << "at position [" << i << "," << j << "]\n";
        }
    }

    // default constructor.
    Matrix B;
    EXPECT_EQ(0, B.row());
    EXPECT_EQ(0, B.col());
}
