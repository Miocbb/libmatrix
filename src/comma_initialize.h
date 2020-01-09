#ifndef _MATRIX_COMMA_INITIALIZE_H_
#define _MATRIX_COMMA_INITIALIZE_H_

#include "matrix.h"

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
            sig_err("Error in MatrixCommaInitializer constructor: trying to initialize a matrix that is not allocated.");
        }
        matrix_(0, 0) = a;
    }

    ~MatrixCommaInitializer()
    {
        if (counter_ < matrix_.size()) {
            std::cout << "Error in matrix with comma initialization: too few elements.";
            std::exit(EXIT_FAILURE);
        } else if (counter_ > matrix_.size()) {
            std::cout << "Error in matrix with comma initialization: too many elements.";
            std::exit(EXIT_FAILURE);
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
            std::cout << "Error in matrix with comma initialization: too many elements.";
            std::exit(EXIT_FAILURE);
        }
        this->matrix_.data()[this->counter_ - 1] = a;
        return *this;
    }
};

}

#endif // _MATRIX_COMMA_INITIALIZE_H_
