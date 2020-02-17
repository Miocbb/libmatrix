/**
 * @file
 *
 * @brief The only header file to include to use the matrix library.
 *
 * @details The purpose of this library is to let the usage of blas/lapack library
 *  for matrice in an easier way.
 *
 * @note Matrix storage format is always the column-wise full storage.
 * @note Matrix index always starts from zero.
 */
#ifndef _MATRIX_INCLUDE_MATRIX_MATRIX_H_
#define _MATRIX_INCLUDE_MATRIX_MATRIX_H_

#include "../src/matrix.h"
#include "../src/matrix_io.h"
#include "../src/comma_initialize.h"
#include "../src/blas.h"
#include "../src/lapack.h"
#include "../src/exception.h"

#endif
