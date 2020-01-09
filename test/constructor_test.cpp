#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

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
    EXPECT_FALSE(A.is_square());
    // check zero initialization.
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < A.col(); j++) {
            EXPECT_DOUBLE_EQ(0.0, A(i, j)) << "Matrix zero initialization not working "
                << "at position [" << i << "," << j << "]\n";
        }
    }

    // default constructor.
    Matrix B;
    EXPECT_TRUE(B.is_square());
    EXPECT_EQ(0, B.row());
    EXPECT_EQ(0, B.col());

    // create a matrix from a vector.
    vector<double> data = {1, 2, 3, 4};
    Matrix C(2, 2, data);
    for (size_t i = 0; i < 4; i++) {
        EXPECT_EQ(C.data()[i], data[i]);
    }
}

TEST(MatrixConstructorTest, comma_initializer)
{
    // matrix A: 1x1.
    Matrix A(1, 1);
    A << 1;
    EXPECT_EQ(A(0,0), 1);

    // matrix B: 2x2;
    Matrix A22(2, 2);
    A22 << 1, 2,
           3, 4;
    for (size_t i = 0; i < 4; i++) {
        EXPECT_EQ(A22.data()[i], i+1);
    }

    //Matrix A00; // run with error
    //A00 << 1;

    //A22 << 1, 2, 3, 4, 5; // run with error.

    //A22 << 1, 2, 3; // run with error.
}
