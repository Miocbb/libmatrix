/**
 * @file matrix.h
 *
 * @brief The interface of matrix library.
 *
 * @details The purpose of this library is to let the usage of blas/lapack library
 *  for matrice in an easier way.
 *
 * @note Matrix storage format is always the column-wise full storage.
 * @note Matrix index always starts from zero.
 */
#ifndef _MATRIX_INCLUDE_MATRIX_MATRIX_H_
#define _MATRIX_INCLUDE_MATRIX_MATRIX_H_

#include "details/matrix.h"
#include "details/matrix_io.h"
#include "details/comma_initialize.h"
#include "details/blas.h"
#include "details/lapack.h"
#include "details/exception.h"

#endif
