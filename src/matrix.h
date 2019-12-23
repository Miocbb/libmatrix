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
#include <random>
#include <cmath>
#include "io.h"
#include "blas_base.h"


namespace matrix {

using std::vector;
using std::size_t;
using std::string;

static std::mt19937 g_rand_generator_mt19937_seed_fixed(1);

class Matrix
{
    private:
    size_t row_;
    size_t col_;
    size_t size_;
    vector<double> data_;

    public:
    /**
     * constructor.
     *
     * the memory for a matrix with input size will
     * be allocated and initialize all matrix elements
     * to be zero.
     *
     * @ param row row number of matrix
     * @ param col column number of matrix
     */
    Matrix(size_t row, size_t col)
    : row_(row), col_(col), data_(row * col, 0.0), size_(row * col) {}

    /**
     * default constructor.
     *
     * creat an empty matrix object.
     */
    Matrix() : row_(0), col_(0), size_(0) {}

    /**
     * access/modify the matrix element by index without
     * bound check.
     *
     * @ param i the row index of the element
     * @ param j the column index of the element
     * @ return referent to the matrix element.
     */
    double & operator()(size_t i, size_t j) {return data_[i * col_ + j];};

    /**
     * access the matrix element by index without bound check.
     *
     * @ param i the row index of the element
     * @ param j the column index of the element
     * @ return const referent to the element value.
     */
    const double & operator()(size_t i, size_t j) const {return data_[i * col_ + j];};

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
    double * data() {return data_.data();}

    /**
     * get the const pointer that points to the begining of
     * matrix data.
     */
    const double * data() const {return data_.data();}

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
     * Make the matrix to be symmetric. The lower triangular
     * matrix data is used by default.
     *
     * @ param use_lower if use the lowertriangular part to
     * symmetrize the matrix. Default it is true.
     * @ return *this the symmetrized matrix itself.
     */
    const Matrix & to_symmetric(bool use_lower = true)
    {
        if (row_ != col_) {
            sig_err("Error: fail to symmetrize the matrix, which isn't square.");
        }
        if (use_lower) {
            for (size_t i = 0; i < row_; i++) {
                for (size_t j = 0; j < i; j++) {
                    (*this)(j, i) = (*this)(i, j);
                }
            }
        } else {
            for (size_t i = 0; i < row_; i++) {
               for (size_t j = 0; j < i; j++) {
                   (*this)(i, j) = (*this)(j, i);
               }
           }
        }
        return *this;
    }

    /**
     * check if the current matrix is symmetric or not based on
     * input accuracy. Default accuracy is 1e-10.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_symmetric(double threshold = 1e-10) const
    {
        threshold = std::fabs(threshold);
        if (row_ != col_) {
            return false;
        }
        const Matrix & T = *this;
        for (size_t i = 0; i < row_; i++) {
            for (size_t j = 0; j < i; j++) {
                if (std::fabs(T(i, j) - T(j, i)) > threshold) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * check if the current matrix is a diagonal matrix or not based on
     * input accuracy. Default accuracy is 1e-10.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_diagonal(double threshold = 1e-10) const
    {
        threshold = std::fabs(threshold);
        if (row_ != col_) {
            return false;
        }
        const Matrix & T = *this;
        for (size_t i = 0; i < row_; i++) {
            for (size_t j = 0; j < i; j++) {
                if (std::fabs(T(i, j)) > threshold || std::fabs(T(j, i)) > threshold) {
                    return false;
                }
            }
        }
        return true;
    }


    /**
     * check if the current matrix is identity or not based on
     * input accuracy. Default accuracy is 1e-10.
     *
     * @ param [in] threshold the threshold of testing float number's equality.
     * @ return bool.
     */
    bool is_identity(double threshold = 1e-10) const
    {
        threshold = std::fabs(threshold);
        if (fabs(threshold) >= 1) {
            sig_err("Error to test matrix identity: too big threshold to make sense.\n");
        }
        if (row_ != col_) {
            return false;
        }
        const Matrix & T = *this;
        for (size_t i = 0; i < row_; i++) {
            if (std::fabs(T(i, i) - 1.0) > threshold) {
                return false;
            }
        }
        if (! this->is_diagonal()) {
            return false;
        }
        return true;
    }

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
    const Matrix & randomize(double a, double b)
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(a, b);
        for (size_t i = 0; i < this->row(); i++) {
            for (size_t j = 0; j < this->col(); j++) {
                (*this)(i, j) = dis(gen);
            }
        }
        return *this;
    }

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
    const Matrix & randomize_seed_fixed(double a, double b)
    {
        std::uniform_real_distribution<> dis(a, b);
        for (size_t i = 0; i < this->row(); i++) {
            for (size_t j = 0; j < this->col(); j++) {
                (*this)(i, j) = dis(g_rand_generator_mt19937_seed_fixed);
            }
        }
        return *this;
    }

    /**
     * Scales current matrix by a constant. A = alpha * A.
     *
     * @ param[in] alpha the scalar coefficient.
     * @ return *this the scaled matrix itself.
     */
    const Matrix & scale(const double alpha)
    {
        int size = size_;
        blas::dscal_(&size, &alpha, data_.data(), blas::ione);
        return *this;
    }

    /**
     * print out the full matrix.
     */
    void show_full() const
    {
        printf("dimension: %zu x %zu, showing in full.\n", row_, col_);
        const size_t numCol = 5;
        size_t k = 0;
        for (size_t i = 0; i < row_; i++) {
            printf(" %5zu:\n", i + 1);
            for (int j = 1; j <= col_; j++) {
                printf(" %15.8e,", data_[k]);
                if (j % numCol == 0 && j != col_) {
                    printf("\n");
                }
                k++;
            }
            printf("\n");
        }
        fflush(stdout);
    }

    /**
     * print out the lower triangular matrix (not the restricted lower part).
     */
    void show_lower() const
    {
        printf("dimension: %zu x %zu, showing the lower triangular parts.\n", row_, col_);
        const size_t numCol = 5;
        for (size_t i = 0; i < row_; i++) {
            printf(" %5zu:\n", i + 1);
            for (int j = 0; j <= i; j++) {
                int ij = i * col_ + j;
                printf(" %15.8e,", data_[ij]);
                if ((j + 1) % numCol == 0 && j != i) {
                    printf("\n");
                }
            }
            printf("\n");
        }
        printf("\n");
        fflush(stdout);
    }

};

}

#endif
