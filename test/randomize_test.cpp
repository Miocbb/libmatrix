#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

/**
 * Test matrix::Matrix::randomize() method.
 *
 * 1. the matrix elements are uniformly distributed in range [a, b).
 * 2. Calling this method each time should give a different random matrix.
 * 3. The seed of the random number generator is not fixed,
 * that is, the random behaviour of this method is not repeatable
 * when run the executable after compilation.
 */

/**
 * Test property #1
 */
TEST(MatrixRandomizeTest, value_range_test)
{
    Matrix A(100, 1000);
    A.randomize(0, 1);
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < A.col(); j++) {
            EXPECT_TRUE((A(i, j) >= 0) && (A(i, j) < 1)) << "Random number is out of range [0, 1) at"
                << " position [" << i << "," << j << "]\n";
        }
    }
}

/**
 * Test property #3.
 */
TEST(MatrixRandomizeTest, repeat_calling_test)
{
    Matrix A(100, 1000);
    A.randomize(0, 1);
    Matrix A_before = A;
    A.randomize(0, 1);
    for (size_t i = 0; i < A.row(); i++) {
        for(size_t j = 0; j < A.col(); j++) {
            EXPECT_NE(A(i, j), A_before(i, j)) << "Warning: random A(i,j) is not changed at i="
                << i << ", j=" << j << ". Check if it is coincidence.\n";
        }
    }
}

TEST(RandomMatrixTest, randoness_eyeball_test)
{
    Matrix A(2, 3);
    std::cout << "Here is a random general matrix in range [0, 1). Eyeball checking. Does it look alright?\n";
    std::cout << "Repeating runing test executable, you should always get a different random matrix.\n";
    A.randomize(0, 1);
    A.show_full();
    A.show_lower();
}
