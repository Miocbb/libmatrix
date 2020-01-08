/**
 * file: matrix.h
 */

#include "matrix.h"
#include "blas_base.h"
#include <random>
#include <cmath>

namespace matrix {

static std::mt19937 g_rand_generator_mt19937_seed_fixed(1);

/**
 * print out the full matrix.
 */
void Matrix::show_full() const
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
void Matrix::show_lower() const
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

/**
 * check if the current matrix is symmetric or not based on
 * input threshold. Default threshold is 1e-10.
 *
 * @ param [in] threshold the threshold of testing float number's equality.
 * @ return bool.
 */
bool Matrix::is_symmetric(double threshold) const
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
 * input threshold. Default threshold is 1e-10.
 *
 * @ param [in] threshold the threshold of testing float number's equality.
 * @ return bool.
 */
bool Matrix::is_diagonal(double threshold) const
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
 * input threshold. Default threshold is 1e-10. The threshold
 * has to be equal or less than 1e-3 to make sense.
 *
 * @ param [in] threshold the threshold of testing float number's equality.
 * @ return bool.
 */
bool Matrix::is_identity(double threshold) const
{
    threshold = std::fabs(threshold);
    if (fabs(threshold) >= 1e-3) {
        sig_err("Error to test matrix identity: threshold is too big to make sense.\n");
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
    if (!this->is_diagonal(threshold)) {
        return false;
    }
    return true;
}

/**
 * check if the current matrix is a zero matrix or not based on
 * input threshold. Default threshold is 1e-10. The threshold has
 * to be equal or less than 1e-3 to make sense.
 *
 * @ param [in] threshold the threshold of testing float number's equality.
 * @ return bool.
 */
bool Matrix::is_zeros(double threshold) const
{
    threshold = std::fabs(threshold);
    if (fabs(threshold) >= 1e-3) {
        sig_err("Error to test zero matrix: threshold is too big to make sense.\n");
    }
    const Matrix & T = *this;
    for (size_t i = 0; i < size_; i++) {
        if (std::fabs(T.data()[i]) > threshold) {
            return false;
        }
    }
    return true;
}

/**
 * check if two matrix is equal by scanning everything based on
 * a threshold. By default, the threshold is 1e-10.
 *
 * @ param A the compared matrix.
 * @ return bool
 */
bool Matrix::is_equal_to(Matrix & A, double threshold) const
{
    threshold = std::fabs(threshold);
    if (size_ != A.size() || row_ != A.row() || col_ != A.col()) {
        return false;
    }
    for (size_t i = 0; i < size_; i++) {
        if (std::fabs(data_[i] - A.data()[i]) > threshold) {
            return false;
        }
    }
    return true;
}

/**
 * calculate matrix trace.
 */
double Matrix::trace() const
{
    if (!this->is_square()) {
        sig_err("Error to get matrix trace: matrix is not square.");
    }
    double rst = 0.0;
    for (size_t i = 0; i < row_; i++) {
        rst += (*this)(i, i);
    }
    return rst;
}

/**
 * Make the matrix to be symmetric. The lower triangular
 * matrix data is used by default.
 *
 * @ param use_lower if use the lowertriangular part to
 * symmetrize the matrix. Default it is true.
 * @ return *this the symmetrized matrix itself.
 */
Matrix & Matrix::to_symmetric(bool use_lower)
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
 * Make the matrix to be random with elements uniformly distributed in range [a, b).
 *
 * THe random number generator is initialized with a non-fixed seed. So the randomness behavior
 * is not repeatable at running time.
 *
 * @ parameter a left range bound.
 * @ parameter b right range bound.
 * @ return *this the random matrix itself.
 */
Matrix & Matrix::randomize(double a, double b)
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
Matrix & Matrix::randomize_seed_fixed(double a, double b)
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
Matrix & Matrix::scale(const double alpha)
{
    int size = size_;
    blas::dscal_(&size, &alpha, data_.data(), blas::ione);
    return *this;
}

/**
 * fill all elements with input number.
 *
 * @ param[in] a the number to be filled.
 * @ return *this the matrix itself.
 */
Matrix & Matrix::fill_all(double a)
{
    #ifdef DOPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < size_; i++) {
        data_[i] = a;
    }
}

}
