/**
 * @file exception.h
 * @brief declaration relate exceptions used in matrix library.
 */
#ifndef _MATRIX_INCLUDE_MATRIX_DETAILS_EXCEPTION_H_
#define _MATRIX_INCLUDE_MATRIX_DETAILS_EXCEPTION_H_

#include "matrix.h"
#include <sstream> // std::stringstream
#include <stdexcept>
#include <string>

namespace matrix {

/**
 * @brief Top-level exception namespace for matrix library.
 */
namespace exception {

using std::string;
using std::stringstream;

/**
 * @brief matrix exception class base.
 */
class MatrixException : public std::runtime_error {
  private:
    void make_message(const char *msg);

  protected:
    /**
     * @brief Exception message.
     */
    stringstream msg_;

  public:
    /**
     * @param [in] msg: Descriptive message for the exception.
     */
    MatrixException(const std::string &msg) : std::runtime_error(msg)
    {
        make_message(msg.c_str());
    }

    /**
     * @brief Get the `const char*` pointer that points to the exception message
     * string.
     */
    const char *what() const noexcept override { return msg_.str().c_str(); }
};

/**
 * @brief Matrix exception class relates to index range errors.
 */
class IndexRangeError : public MatrixException {
  public:
    /**
     * @param [in] msg: Descriptive message for the exception.
     */
    IndexRangeError(const string &msg)
        : MatrixException("Matrix index range error")
    {
        msg_ << "Description:" << msg << std::endl;
    }
};

/**
 * @brief Matrix exception class relates to matrix dimemsion errors.
 */
class DimensionError : public MatrixException {
  public:
    /**
     * @brief Create an exception that the matrix dimension between two matrices
     * is not matched.
     * @param [in] A1: The first matrix.
     * @param [in] A2: The first matrix.
     * @param [in] msg: Detailed descriptive message for the exception.
     */
    DimensionError(const Matrix &A1, const Matrix &A2, const string &msg)
        : MatrixException("Two matrices dimension not matched.")
    {
        msg_ << "Description: " << msg << std::endl;
        msg_ << "Details: "
             << "Matrix 1 dimension: [" << A1.row() << ", " << A1.col() << "]"
             << std::endl;
        msg_ << "         "
             << "Matrix 2 dimension: [" << A2.row() << ", " << A2.col() << "]"
             << std::endl;
    }

    /**
     * @brief Create an exception for the general unmatched dimension error with
     * size information.
     * @param [in] expected: The expected dimension.
     * @param [in] actual: The actual dimension.
     * @param [in] msg: Detailed descriptive message for the exception.
     */
    DimensionError(size_t expected, size_t actual, const string &msg)
        : MatrixException("Dimension does not match with the expectation.")
    {
        msg_ << "Description: " << msg << std::endl;
        msg_ << "Details: "
             << "actual dimension is " << actual
             << ", while the expected one should be " << expected << "."
             << std::endl;
    }

    /**
     * @brief Create an exception for the general dimension error.
     * @param [in] msg: Detailed descriptive message for the exception.
     */
    DimensionError(const string &msg) : MatrixException("Dimension error.")
    {
        msg_ << "Description: " << msg << std::endl;
    }
};

/**
 * @brief Matrix exception class relates to matrix operation errors.
 */
class MatrixOperationError : public MatrixException {
  public:
    /**
     * @brief Create an exception for the matrix operation error with operation
     * function name.
     * @param [in] op_func: The name of the function that perform the matrix
     * operation.
     * @param [in] msg: Detailed descriptive message for the exception.
     */
    MatrixOperationError(const string &op_func, const string &msg)
        : MatrixException("Matrix operation error")
    {
        msg_ << "Matrix operation name: " << op_func << std::endl;
        msg_ << "Details: " << msg << std::endl;
    }
};

/**
 * @brief Matrix exception class relates to matrix I/O errors.
 */
class MatrixIOException : public MatrixException {
  public:
    /**
     * @brief Create an exception for the matrix I/O error with file name.
     * @param [in] file: The name of file relates to the I/O error.
     * @param [in] msg: Detailed descriptive message for the exception.
     */
    MatrixIOException(const string &file, const string &msg)
        : MatrixException("Matrix I/O error")
    {
        msg_ << "I/O file name: " << file << std::endl;
        msg_ << "Details: " << msg << std::endl;
    }
};

} // namespace exception
} // namespace matrix

#endif // _MATRIX_SRC_EXCEPTION_H_
