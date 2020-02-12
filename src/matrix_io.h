#ifndef _MATRIX_MATRIX_IO_H_
#define _MATRIX_MATRIX_IO_H_

#include "matrix.h"
#include <memory>

namespace matrix {

/**
 * write a number of matrix into binary file in order.
 */
void write_matrices_to_binary(vector<std::shared_ptr<const Matrix>> &Mat, const char *fname);

/**
 * read a number of matrices from a binary file.
 */
void read_matrices_from_binary(vector<std::shared_ptr<Matrix>> &Mat, const char *fname);

}

#endif  // _MATRIX_MATRIX_IO_H_
