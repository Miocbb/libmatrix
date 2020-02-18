#include <gtest/gtest.h>
#include <matrix/matrix.h>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

using matrix::Matrix;
TEST(MatrixWriteReadTxtTest, write_read_pair_test)
{
    std::vector<std::shared_ptr<Matrix>> mat;
    mat.push_back(std::make_shared<Matrix>(2, 2));
    // one matrix.
    *mat[0] = {1, 2, 3, 4};
    std::string file_path = realpath(__FILE__, NULL);
    std::string txt_path =
        file_path.substr(0, file_path.rfind("/")) + "/test.csv.tem";
    matrix::write_matrices_to_txt(mat, txt_path);
    auto mat_read = matrix::read_matrices_from_txt(txt_path);
    for (size_t i = 0; i < mat_read.size(); i++) {
        auto status = mat_read[i]->is_equal_to(*mat[i]);
        if (status == false) {
            std::cout << "Matrix read: i = " << i + 1 << std::endl;
            mat_read[i]->show_full();
        }
        EXPECT_TRUE(status);
    }

    // three matrices.
    mat.push_back(std::make_shared<Matrix>(3, 3));
    mat.push_back(std::make_shared<Matrix>(2, 3));
    *mat[1] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    *mat[2] = {1, 2, 3, 1, 2, 3};
    matrix::write_matrices_to_txt(mat, txt_path);
    mat_read = matrix::read_matrices_from_txt(txt_path);
    for (size_t i = 0; i < mat_read.size(); i++) {
        auto status = mat_read[i]->is_equal_to(*mat[i]);
        if (status == false) {
            std::cout << "Matrix read: i = " << i + 1 << std::endl;
            mat_read[i]->show_full();
        }
        EXPECT_TRUE(status);
    }
}
