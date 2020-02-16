#ifndef _MATRIX_SRC_MATRIX_IO_H_
#define _MATRIX_SRC_MATRIX_IO_H_

#include <memory>
#include "matrix.h"

namespace matrix {

/**
 * write a number of matrix into binary file in order.
 */
void write_matrices_to_binary(vector<std::shared_ptr<const Matrix>> &Mat, const char *fname);

/**
 * read a number of matrices from a binary file.
 */
void read_matrices_from_binary(vector<std::shared_ptr<Matrix>> &Mat, const char *fname);

/**
 * read matrix from binary file and create a vector of matrix.
 */
std::vector<std::shared_ptr<Matrix>> read_matrices_from_binary(const char *fname);
std::vector<std::shared_ptr<Matrix>> read_matrices_from_binary(string &fname);

}   // namespace matrix

#endif  // _MATRIX_SRC_MATRIX_IO_H_
