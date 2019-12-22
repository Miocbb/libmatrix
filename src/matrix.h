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
#include "io.h"


namespace matrix {

using std::vector;
using std::size_t;
using std::string;

class Matrix
{
    private:
    size_t row_;
    size_t col_;
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
    : row_(row), col_(col), data_(row * col, 0.0) {}

    /**
     * default constructor.
     *
     * creat an empty matrix object.
     */
    Matrix() : row_(0), col_(0) {}

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
     * make the matrix to be random with elements uniformly distributed in range [a, b).
     *
     * @ parameter a left range bound.
     * @ parameter b right range bound.
     * @ return *this the random matrix itself.
     *
     * 1. the matrix elements are uniformly distributed in range [a, b).
     * 2. Calling this method each time should give a different random matrix.
     * 3. The seed of the random number generator is not fixed,
     * that is, the random behaviour of this method is not repeatable
     * when run the executable after compilation.
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
