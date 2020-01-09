#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

TEST(MatrixConstructorTest, index_only_constructor)
{
    Matrix A(2, 3);
    // check zero initialization.
    for (size_t i = 0; i < A.row(); i++) {
        for (size_t j = 0; j < A.col(); j++) {
            EXPECT_DOUBLE_EQ(0.0, A(i, j)) << "Matrix zero initialization not working "
                << "at position [" << i << "," << j << "]\n";
        }
    }
}

TEST(MatrixConstructorTest, construct_from_vector)
{
    vector<double> data = {1, 2, 3, 4, 5, 6};
    // deep copy from a vector.
    Matrix A(2, 3, data, matrix::Matrix::CopyType::kDeepCopy);
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_DOUBLE_EQ(A.data()[i], data[i]);
    }
    EXPECT_TRUE(!A.is_data_stored_outside());
    // shallow copy from a vector.
    Matrix B(2, 3, data, matrix::Matrix::CopyType::kShallowCopy);
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_DOUBLE_EQ(B.data()[i], data[i]);
    }
    EXPECT_TRUE(B.is_data_stored_outside());

    //// size matching checking
    //Matrix C(2, 2, data, matrix::Matrix::CopyType::kDeepCopy); // run with error
    //Matrix D(2, 4, data, matrix::Matrix::CopyType::kDeepCopy); // run with error
}

TEST(MatrixConstructorTest, construct_from_pointer)
{
    vector<double> data = {1, 2, 3, 4, 5, 6};
    // deep copy from a pointer.
    Matrix A(2, 3, data.data(), matrix::Matrix::CopyType::kDeepCopy);
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_DOUBLE_EQ(A.data()[i], data[i]);
    }
    EXPECT_TRUE(!A.is_data_stored_outside());
    // shallow copy from a pointer.
    Matrix B(2, 3, data.data(), matrix::Matrix::CopyType::kShallowCopy);
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_DOUBLE_EQ(B.data()[i], data[i]);
    }
    EXPECT_TRUE(B.is_data_stored_outside());

    // smaller size, the first 4 elements should match with the data array.
    Matrix C(2, 2, data.data(), matrix::Matrix::CopyType::kShallowCopy);
    for (size_t i = 0; i < data.size(); i++) {
        EXPECT_DOUBLE_EQ(C.data()[i], data[i]);
    }
    EXPECT_TRUE(C.is_data_stored_outside());
}

TEST(MatrixConstructorTest, default_constructor)
{
    Matrix A;
    EXPECT_EQ(A.row(), 0.0);
    EXPECT_EQ(A.col(), 0.0);
    EXPECT_EQ(A.size(), 0.0);
    EXPECT_TRUE(A.data() == nullptr);
}

TEST(MatrixConstructorTest, copy_constructor)
{
    Matrix A(2, 2);
    Matrix B(A);
    EXPECT_TRUE(B.is_equal_to(A));
    EXPECT_TRUE(B.data() != A.data());  // make sure it is deep copy.
    EXPECT_TRUE(!B.is_data_stored_outside());
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
