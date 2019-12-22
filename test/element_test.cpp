#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <vector>

using matrix::Matrix;
using std::vector;

/**
 * test operator() and .at() methods to access/modify matrix element.
 */
TEST(ElementAccessTest, get_assign_element_test)
{
    vector<double> A23_data = {1, 2, 3, 4, 5, 6};
    Matrix A(2, 3);
    size_t n = 0;
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 3; j++) {
            A(i, j) = A23_data[n];
            n++;
        }
    }
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(A23_data[i], A.data()[i]) << "Error in matrix element accessor (), postion: "
            << i << "\n";
    }

    n = 0;
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 3; j++) {
            A.at(i, j) = A23_data[n];
            n++;
        }
    }
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(A23_data[i], A.data()[i]) << "Error in matrix element accessor .at(), postion: "
            << i << "\n";
    }
}
