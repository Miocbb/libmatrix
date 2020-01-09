/**
 * file: matrix.h
 *
 * matrix class decleration.
 */

#ifndef _MATRIX_MATRIX_H_
#define _MATRIX_MATRIX_H_

#include <vector>
#include <cstddef>
#include <string>
#include <iostream>
#include "io.h"
#include "blas_base.h"


namespace matrix {

using std::vector;
using std::size_t;
using std::string;

/**
 * Helper class used to do comma initialization for matrix object like Eigen library.
 */
class MatrixCommaInitializer;

class Matrix
{
    private:
    size_t row_;
    size_t col_;
    size_t size_;
    /* The data of matrix is either stored in vector `data_` or outside of the
     * object that is pointed by the double pointer `data_ptr_`. */
    vector<double> data_vec_;
    double * data_ptr_;     /* This pointer will always point to the head of the matrix data.
                             * No memory allocation is assigned for this pointer. */

    public:

    enum CopyType {
        kShallowCopy,
        kDeepCopy,
    };

    /**
     * construct a matrix with given size.
     *
     * the memory for a matrix with input size will
     * be allocated and initialize all matrix elements
     * to be zero.
     *
     * matrix data is stored in inside of this object.
     *
     * @ param row row number of matrix
     * @ param col column number of matrix
     */
    Matrix(size_t row, size_t col)
    : row_(row), col_(col), data_vec_(row * col, 0.0), size_(row * col),
    data_ptr_{data_vec_.data()} {}

    /**
     * construct a matrix from a std::vector<double>.
     *
     * The data size of the vector will be check to match the matrix size.
     *
     * This constructor is not const friendly. If you want to do a shallow copy
     * of a const vector to the matrix, you can do declare the created Matrix to
     * be const and typecast the const vector to be non-const vector, that is,
     * ```
     * const vector<double> data;
     * const Matrix A(row, col, const_cast <vector<double>>(data), kShallowCopy);
     * ```
     * This can preserve the created matrix not change the original vector data.
     *
     * @ param[in] data: the input matrix data.
     * @ param[in] copy_type: If copy_type is deep copy, all the data
     *  will be copied into matrix. If the copy type is shallow copy, the pointer to the
     *  data memory is assigned to the matrix and the matrix gains the access to the data.
     */
    Matrix(size_t row, size_t col, vector<double> & inp_data, CopyType copy_type)
    : row_(row), col_(col), size_(row * col), data_vec_(0)
    {
        if (size_ != inp_data.size()) {
            std::cout << "Error in creating a matrix from a vector object: unmatched size.";
            std::exit(EXIT_FAILURE);
        }
        if (copy_type == kDeepCopy) {
            data_vec_ = inp_data;
            data_ptr_ = data_vec_.data();
        } else if (copy_type == kShallowCopy) {
            data_ptr_ = inp_data.data();
        } else {
            sig_err("Error in Matrix constructor: unknown copy type.");
        }
    }

    /**
     * construct a matrix from a double array pointer.
     *
     * The data size of the array will NOT be check to match the matrix size.
     * Use this constructor carefully.
     *
     * This constructor is not const friendly. If you want to do a shallow copy
     * of a const double array to the matrix, you can do declare the created Matrix to
     * be const and typecast the const array to be non-const array, that is,
     * ```
     * const double * data;
     * const Matrix A(row, col, const_cast <double *>(data), kShallowCopy);
     * ```
     * This can preserve the created matrix not change the original array data.
     *
     * @ param[in] data: the input matrix data.
     * @ param[in] copy_type: If copy_type is deep copy, all the data
     *  will be copied into matrix. If the copy type is shallow copy, the pointer to the
     *  data memory is assigned to the matrix and the matrix gains the access to the data.
     */
    Matrix(size_t row, size_t col, double *inp_data_ptr, CopyType copy_type)
    : row_(row), col_(col), size_(row * col), data_vec_(0)
    {
        if (copy_type == kDeepCopy) {
            data_vec_.resize(size_);
            data_ptr_ = data_vec_.data();
            int dim = size_;
            blas::dcopy_(&dim, inp_data_ptr, blas::ione, data_ptr_, blas::ione);
        } else if (copy_type == kShallowCopy) {
            data_ptr_ = inp_data_ptr;
        } else {
            sig_err("Error in Matrix constructor: unknown copy type.");
        }
    }

    /**
     * default constructor.
     *
     * creat an empty matrix object.
     */
    Matrix() : row_(0), col_(0), size_(0), data_vec_(0), data_ptr_{nullptr} {}

    /**
     * Copy constructor.
     *
     * Always do a deep copy.
     */
    Matrix(const Matrix & other) : row_{other.row()}, col_{other.col()}, size_{other.size()},
    data_vec_(size_), data_ptr_{data_vec_.data()}
    {
        int dim = size_;
        blas::dcopy_(&dim, other.data(), blas::ione, data_ptr_, blas::ione);
    }

    /**
     * Copy assignment operator overloading.
     *
     * Always make a deep copy of the matrix and assign it to the destination.
     */
    Matrix & operator = (const Matrix & other)
    {
        row_ = other.row();
        col_ = other.col();
        size_ = other.size();
        data_vec_.resize(size_);
        data_ptr_ = data_vec_.data();
        int dim = size_;
        blas::dcopy_(&dim, other.data(), blas::ione, data_ptr_, blas::ione);
        return *this;
    }

    /**
     * Initialize matrix object with the help of `<<` and `,` operator in an easy way.
     * That is, `A << 1, 2, 3, 4;`. The matrix dimension of A has to be specified
     * in advance and the number of elements has to be equal to the matrix size,
     * otherwise, it will abort with error.
     */
    MatrixCommaInitializer operator <<(double a);

    /**
     * access/modify the matrix element by index without
     * bound check.
     *
     * @ param i the row index of the element
     * @ param j the column index of the element
     * @ return referent to the matrix element.
     */
    double & operator()(size_t i, size_t j) {return data_ptr_[i * col_ + j];};

    /**
     * access the matrix element by index without bound check.
     *
     * @ param i the row index of the element
     * @ param j the column index of the element
     * @ return const referent to the element value.
     */
    const double & operator()(size_t i, size_t j) const {return data_ptr_[i * col_ + j];};

    /**
     * access/modify the matrix element by index with bound check.
     *
     * @ param i the row index of the element
     * @ param j the column index of the element
     * @ return referent to the matrix element.
     */
    double & at(size_t i, size_t j)
    {
        if (i >= row_ || j >= col_) {
            sig_err("Error: matrix index is out of bound.\n");
        }
        return (*this)(i, j);
    }

    /**
     * access the matrix element by index with bound check.
     *
     * @ param i the row index of the element
     * @ param j the column index of the element
     * @ return referent to the matrix element.
     */
    const double & at(size_t i, size_t j) const
    {
        if (i >= row_ || j >= col_) {
            sig_err("Error: matrix index is out of bound.\n");
        }
        return (*this)(i, j);
    }

    /**
     * get the pointer that points to the begining of
     * matrix data.
     */
    double * data() {return data_ptr_;}

    /**
     * get the const pointer that points to the begining of
     * matrix data.
     */
    const double * data() const {return data_ptr_;}

    /**
     * get matrix row numbers.
     */
    const size_t & row() const {return row_;}

    /**
     * get matrix column numbers.
     */
    const size_t & col() const {return col_;}

    /**
     * get matrix size, that is, how many elements the matrix has.
     */
    const size_t & size() const {return size_;}

    /**
     * check if the data is stored in vector or outside of the object.
     *
     * special case:
     * if the matrix is an empty matrix, it will return false.
     */
    bool is_data_stored_outside() const
    {
        if (this->size() == 0) {
            return false;
        } else {
            return data_ptr_ != data_vec_.data();
        }
    }

    /**
     * check if the matrix is square or not.
     */
    bool is_square() const {return (row_ == col_);}

    /**
     * check if the current matrix is symmetric or not based on
     * input threshold. Default threshold is 1e-10.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_symmetric(double threshold = 1e-10) const;

    /**
     * check if the current matrix is a diagonal matrix or not based on
     * input threshold. Default threshold is 1e-10.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_diagonal(double threshold = 1e-10) const;

    /**
     * check if the current matrix is identity or not based on
     * input threshold. Default threshold is 1e-10. The threshold
     * has to be equal or less than 1e-3 to make sense.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_identity(double threshold = 1e-10) const;

    /**
     * check if the current matrix is a zero matrix or not based on
     * input threshold. Default threshold is 1e-10. The threshold has
     * to be equal or less than 1e-3 to make sense.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_zeros(double threshold = 1e-10) const;

    /**
     * check if two matrix is equal by scanning everything based on
     * a threshold. By default, the threshold is 1e-10.
     *
     * @ param A the compared matrix.
     * @ return bool
     */
    bool is_equal_to(Matrix & A, double threshold = 1e-10) const;

    /**
     * print out the full matrix.
     */
    void show_full() const;

    /**
     * print out the lower triangular matrix (not the restricted lower part).
     */
    void show_lower() const;

    /**
     * calculate matrix trace.
     */
    double trace() const;

    /**
     * Make the matrix to be symmetric.
     *
     * @ param[in] uplo "U" using upper triangular part; "L" using lower triangular part.
     * @ return *this the symmetrized matrix itself.
     */
    Matrix & to_symmetric(const string & uplo);

    /**
     * Make the matrix to be random with elements uniformly distributed in range [a, b).
     *
     * THe random number generator is initialized with a non-fixed seed. So the randomness behavior
     * is not repeatable at running time.
     *
     * @ parameter a left range bound.
     * @ parameter b right range bound.
     * @ return *this the random matrix itself.
     */
    Matrix & randomize(double a, double b);

    /**
     * Make the matrix to be random with elements uniformly distributed in range [a, b).
     *
     * THe random number generator is initialized with a fixed seed. So the randomness behavior
     * is repeatable at running time.
     *
     * @ parameter a left range bound.
     * @ parameter b right range bound.
     * @ return *this the random matrix itself.
     */
    Matrix & randomize_seed_fixed(double a, double b);

    /**
     * Scales current matrix by a constant. A = alpha * A.
     *
     * @ param[in] alpha the scalar coefficient.
     * @ return *this the scaled matrix itself.
     */
    Matrix & scale(const double alpha);

    /**
     * fill all elements with input number.
     *
     * @ param[in] a the number to be filled.
     * @ return *this the matrix itself.
     */
    Matrix & fill_all(double a);

};

}

#endif
