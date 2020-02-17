/**
 * @file
 * @brief declarations relate to matrix I/O.
 */
#ifndef _MATRIX_SRC_MATRIX_IO_H_
#define _MATRIX_SRC_MATRIX_IO_H_

#include "matrix.h"
#include <memory>

namespace matrix {

/**
 * @brief Write a number of matrix into binary file in order.
 * @param [in] Mat: a vector matrices to be written.
 * @param [in] fname: the binary file name (relative/absolute path).
 */
void write_matrices_to_binary(vector<std::shared_ptr<const Matrix>> &Mat,
                              const char *fname);

/**
 * @brief Read a number of matrix into binary file in order.
 * @param [in] Mat: a vector matrices to be written.
 * @param [in] fname: the binary file name (relative/absolute path).
 */
void read_matrices_from_binary(vector<std::shared_ptr<Matrix>> &Mat,
                               const char *fname);

/**
 * @brief Read matrix/matrices from binary file and create a vector of matrices.
 * @param [in] fname: binary file name.
 * @return std::vector<std::shared_ptr<Matrix>>: a vector of matrices.
 */
std::vector<std::shared_ptr<Matrix>>
read_matrices_from_binary(const char *fname);

/**
 * @brief Read matrix/matrices from binary file and create a vector of matrices.
 * @param [in] fname: binary file name.
 * @return std::vector<std::shared_ptr<Matrix>>: a vector of matrices.
 */
std::vector<std::shared_ptr<Matrix>> read_matrices_from_binary(string &fname);

} // namespace matrix

#endif // _MATRIX_SRC_MATRIX_IO_H_
