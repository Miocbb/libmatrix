#ifndef _MATRIX_EXCEPTION_H_
#define _MATRIX_EXCEPTION_H_

#include <stdexcept>
#include <string>
#include <sstream>  // std::stringstream
#include "matrix.h"

namespace matrix {

namespace exception {

using std::string;
using std::stringstream;

class MatrixException : public std::runtime_error
{
private:
    void make_message(const char* msg);
protected:
    stringstream msg_;
public:
    MatrixException(const std::string& msg) : std::runtime_error(msg) {make_message(msg.c_str());}

    const char* what() const noexcept override {return msg_.str().c_str();}
};

class IndexRangeError : public MatrixException
{
public:
    IndexRangeError(const string& msg)
        : MatrixException("Matrix index range error")
    {
        msg_ << "Description:" << msg << std::endl;
    }
};

class DimensionError : public MatrixException
{
public:
    DimensionError(const Matrix &A1, const Matrix &A2, const string& msg)
        : MatrixException("Two matrices dimension not matched.")
    {
        msg_ << "Description: " << msg << std::endl;
        msg_ << "Details: " << "Matrix 1 dimension: [" << A1.row() << ", " << A1.col() << "]" << std::endl;
        msg_ << "         " << "Matrix 2 dimension: [" << A2.row() << ", " << A2.col() << "]" << std::endl;
    }

    DimensionError(size_t expect, size_t input, const string& msg)
        : MatrixException("Dimension does not match with the expectation.")
    {
        msg_ << "Description: " << msg << std::endl;
        msg_ << "Details: " << "Dimension is " << input
             << ", while the expected one should be " << expect << "." << std::endl;
    }

    DimensionError(const string& msg) : MatrixException("Dimension error.")
    {
        msg_ << "Description: " << msg << std::endl;
    }
};

class MatrixOperationError : public MatrixException
{
public:
    MatrixOperationError(const string& op_func, const string& msg)
        : MatrixException("Matrix operation error")
    {
        msg_ << "Matrix operation name: " << op_func << std::endl;
        msg_ << "Details: " << msg << std::endl;
    }
};

class MatrixIOException : public MatrixException
{
public:
    MatrixIOException(const string& file, const string& msg)
        : MatrixException("Matrix I/O error")
    {
        msg_ << "I/O file name: " << file << std::endl;
        msg_ << "Details: " << msg << std::endl;
    }
};


}

}

#endif // _MATRIX_EXCEPTION_H_
