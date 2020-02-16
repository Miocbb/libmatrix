#ifndef _MATRIX_SRC_COMMA_INITIALIZE_H_
#define _MATRIX_SRC_COMMA_INITIALIZE_H_

#include <string>
#include "matrix.h"
#include "exception.h"

namespace matrix {
/**
 * Helper class used to do comma initialization for matrix object like Eigen library.
 */
class MatrixCommaInitializer
{
    private:
    size_t counter_;
    Matrix & matrix_;

    public:
    MatrixCommaInitializer(Matrix & A, double a) : counter_{1}, matrix_{A}
    {
        if (matrix_.size() == 0) {
            throw exception::MatrixException("Error in MatrixCommaInitializer constructor: trying to initialize a matrix that is not allocated.");
        }
        matrix_(0, 0) = a;
    }

    ~MatrixCommaInitializer()
    {
        if (counter_ < matrix_.size()) {
            string msg{"Error in `matrix::Matrix` with comma initialization: too few elements."};
            throw exception::DimensionError(matrix_.size(), counter_, msg);
        } else if (counter_ > matrix_.size()) {
            string msg{"Error in `matrix::Matrix` with comma initialization: too many elements."};
            throw exception::DimensionError(matrix_.size(), counter_, msg);
        }
    }

    /**
     * overloading ',' operator to insert data into matrix. The number of inserted data element
     * is tracked. When too many data is inserted, it will abort with error.
     */
    MatrixCommaInitializer & operator , (double a)
    {
        this->counter_++;
        if (this->counter_ > this->matrix_.size()) {
            string msg{"Error in `matrix::Matrix` with comma initialization: too many elements."};
            throw exception::DimensionError(matrix_.size(), counter_, msg);
        }
        this->matrix_.data()[this->counter_ - 1] = a;
        return *this;
    }
};

}   // namespace matrix

#endif // _MATRIX_SRC_COMMA_INITIALIZE_H_
