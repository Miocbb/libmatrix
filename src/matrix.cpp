/**
 * @file
 * @brief Definition matrix class.
 */

#include "matrix.h"

#include "blas_base.h"
#include "comma_initialize.h"
#include "exception.h"
#include <cmath>
#include <random>

namespace matrix {

static std::mt19937 g_rand_generator_mt19937_seed_fixed(1);

Matrix::Matrix(size_t row, size_t col, vector<double> &inp_data,
               CopyType copy_type)
    : row_(row), col_(col), size_(row * col), data_vec_(0)
{
    if (size_ != inp_data.size()) {
        throw exception::DimensionError(
            size_, inp_data.size(),
            "Fail to create a `matrix::Matrix` from a `std::vector`.");
    }
    if (copy_type == kDeepCopy) {
        data_vec_ = inp_data;
        data_ptr_ = data_vec_.data();
    } else if (copy_type == kShallowCopy) {
        data_ptr_ = inp_data.data();
    } else {
        throw exception::MatrixException(
            "Fail to create a `matrix::Matrix` object: unknown copy type.");
    }
}

Matrix::Matrix(size_t row, size_t col, double *inp_data_ptr, CopyType copy_type)
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
        throw exception::MatrixException(
            "Fail to create a `matrix::Matrix` object: unknown copy type.");
    }
}

/**
 * @note This will always do the deep copy.
 */
Matrix::Matrix(const Matrix &other)
    : row_{other.row()}, col_{other.col()}, size_{other.size()},
      data_vec_(size_), data_ptr_{data_vec_.data()}
{
    int dim = size_;
    blas::dcopy_(&dim, other.data(), blas::ione, data_ptr_, blas::ione);
}

/**
 * @note This will always make a deep copy of the matrix and assign it to the
 * destination.
 */
Matrix &Matrix::operator=(const Matrix &other)
{
    if (&other == this) {
        return *this;
    }
    row_ = other.row();
    col_ = other.col();
    size_ = other.size();
    data_vec_.resize(size_);
    data_ptr_ = data_vec_.data();
    int dim = size_;
    blas::dcopy_(&dim, other.data(), blas::ione, data_ptr_, blas::ione);
    return *this;
}

const Matrix &Matrix::operator=(std::initializer_list<double> init_list)
{
    if (init_list.size() != this->size()) {
        string msg{"Fail to initialize `matrix::Matrix` with initializer_list "
                   "syntax. Unmatched size."};
        throw exception::DimensionError(this->size(), init_list.size(), msg);
    }
    size_t i = 0;
    for (auto p = init_list.begin(); p != init_list.end(); p++) {
        this->data()[i] = *p;
        i++;
    }
}

const double &Matrix::at(size_t i, size_t j) const
{
    if (i >= row_ || j >= col_) {
        std::stringstream msg;
        msg << "Index is out of range at position (" << i << ", " << j << ".";
        throw exception::IndexRangeError(msg.str());
    }
    return (*this)(i, j);
}

MatrixCommaInitializer Matrix::operator<<(double a)
{
    return MatrixCommaInitializer(*this, a);
}

void Matrix::show_full(size_t elements_per_line) const
{
    printf("dimension: %zu x %zu, showing in full.\n", row_, col_);
    size_t k = 0;
    for (size_t i = 0; i < row_; i++) {
        printf(" %5zu:\n", i + 1);
        for (int j = 1; j <= col_; j++) {
            printf(" %15.8e,", this->data()[k]);
            if (j % elements_per_line == 0 && j != col_) {
                printf("\n");
            }
            k++;
        }
        printf("\n");
    }
    fflush(stdout);
}

void Matrix::show_lower(size_t elements_per_line) const
{
    printf("dimension: %zu x %zu, showing the lower triangular parts.\n", row_,
           col_);
    for (size_t i = 0; i < row_; i++) {
        printf(" %5zu:\n", i + 1);
        for (int j = 0; j <= i; j++) {
            int ij = i * col_ + j;
            printf(" %15.8e,", this->data()[ij]);
            if ((j + 1) % elements_per_line == 0 && j != i) {
                printf("\n");
            }
        }
        printf("\n");
    }
    printf("\n");
    fflush(stdout);
}

bool Matrix::is_symmetric(double threshold) const
{
    threshold = std::fabs(threshold);
    if (row_ != col_) {
        return false;
    }
    const Matrix &T = *this;
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
 * @details It returns true only when all |A(i, j)| < \p threshold,
 *  for any i not equal to j.
 */
bool Matrix::is_diagonal(double threshold) const
{
    threshold = std::fabs(threshold);
    if (row_ != col_) {
        return false;
    }
    const Matrix &T = *this;
    for (size_t i = 0; i < row_; i++) {
        for (size_t j = 0; j < i; j++) {
            if (std::fabs(T(i, j)) > threshold ||
                std::fabs(T(j, i)) > threshold) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::is_identity(double threshold) const
{
    threshold = std::fabs(threshold);
    if (row_ != col_) {
        return false;
    }
    const Matrix &T = *this;
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

bool Matrix::is_zeros(double threshold) const
{
    threshold = std::fabs(threshold);
    const Matrix &T = *this;
    for (size_t i = 0; i < size_; i++) {
        if (std::fabs(T.data()[i]) > threshold) {
            return false;
        }
    }
    return true;
}

bool Matrix::is_equal_to(const Matrix &A, double threshold) const
{
    threshold = std::fabs(threshold);
    if (size_ != A.size() || row_ != A.row() || col_ != A.col()) {
        return false;
    }
    for (size_t i = 0; i < size_; i++) {
        if (std::fabs(this->data()[i] - A.data()[i]) > threshold) {
            return false;
        }
    }
    return true;
}

bool Matrix::is_same_dimension_to(const Matrix &other) const
{
    return ((this->row() == other.row()) && (this->col() == other.col()));
}

double Matrix::trace() const
{
    if (!this->is_square()) {
        throw exception::DimensionError(
            "Cannot get trace of a matrix that is not squared.");
    }
    double rst = 0.0;
    for (size_t i = 0; i < row_; i++) {
        rst += (*this)(i, i);
    }
    return rst;
}

Matrix &Matrix::to_symmetric(const string &uplo)
{
    if (row_ != col_) {
        throw exception::DimensionError(
            "Cannot symmetrize a matrix that is not squared.");
    }
    if (uplo == "U") {
        for (size_t i = 0; i < row_; i++) {
            for (size_t j = 0; j < i; j++) {
                (*this)(i, j) = (*this)(j, i);
            }
        }
    } else if (uplo == "L") {
        for (size_t i = 0; i < row_; i++) {
            for (size_t j = 0; j < i; j++) {
                (*this)(j, i) = (*this)(i, j);
            }
        }
    }
    return *this;
}

Matrix &Matrix::randomize(double a, double b)
{
    std::random_device
        rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < this->row(); i++) {
        for (size_t j = 0; j < this->col(); j++) {
            (*this)(i, j) = dis(gen);
        }
    }
    return *this;
}

Matrix &Matrix::randomize_seed_fixed(double a, double b)
{
    std::uniform_real_distribution<> dis(a, b);
    for (size_t i = 0; i < this->row(); i++) {
        for (size_t j = 0; j < this->col(); j++) {
            (*this)(i, j) = dis(g_rand_generator_mt19937_seed_fixed);
        }
    }
    return *this;
}

Matrix &Matrix::scale(const double alpha)
{
    int size = size_;
    blas::dscal_(&size, &alpha, this->data(), blas::ione);
    return *this;
}

Matrix &Matrix::fill_all(double a)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < size_; i++) {
        this->data()[i] = a;
    }
}

Matrix &Matrix::set_identity()
{
    if (!this->is_square()) {
        throw exception::DimensionError(
            "Cannot make a non-square matrix to be identity.");
    }
    this->fill_all(0.0);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < row_; i++) {
        (*this)(i, i) = 1.0;
    }
}

} // namespace matrix
