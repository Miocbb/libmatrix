/**
 * file: matrix.h
 */

#include "matrix.h"

#include <random>
#include <cmath>
#include "blas_base.h"
#include "comma_initialize.h"
#include "exception.h"

namespace matrix {

static std::mt19937 g_rand_generator_mt19937_seed_fixed(1);

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
Matrix::Matrix(size_t row, size_t col, vector<double> & inp_data, CopyType copy_type)
: row_(row), col_(col), size_(row * col), data_vec_(0)
{
    if (size_ != inp_data.size()) {
        throw exception::DimensionError(size_, inp_data.size(), "Fail to create a `matrix::Matrix` from a `std::vector`.");
    }
    if (copy_type == kDeepCopy) {
        data_vec_ = inp_data;
        data_ptr_ = data_vec_.data();
    } else if (copy_type == kShallowCopy) {
        data_ptr_ = inp_data.data();
    } else {
        throw exception::MatrixException("Fail to create a `matrix::Matrix` object: unknown copy type.");
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
        throw exception::MatrixException("Fail to create a `matrix::Matrix` object: unknown copy type.");
    }
}

/**
 * Copy constructor.
 *
 * Always do a deep copy.
 */
Matrix::Matrix(const Matrix & other)
    : row_{other.row()}, col_{other.col()}, size_{other.size()},
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
Matrix& Matrix::operator = (const Matrix & other)
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

/**
 * Copy assignment operator overloading: enable an easy way to do matrix
 * element initialization with std::initializer_list.
 *
 * The matrix dimension has to be declared in advance. The input list size
 * will be checked for the initialization.
 * e.g.
 * Matrix A(2, 2);
 * A = {1, 2,
 *      3, 4};
 *
 * @ param[in] init_list: the data contained in the initializer_list.
 * @ return: the const reference to the matrix itself.
 */
const Matrix & Matrix::operator = (std::initializer_list<double> init_list)
{
    if (init_list.size() != this->size()) {
        string msg {"Fail to initialize `matrix::Matrix` with initializer_list syntax. Unmatched size."};
        throw exception::DimensionError(this->size(), init_list.size(), msg);
    }
    size_t i = 0;
    for (auto p = init_list.begin(); p != init_list.end(); p++) {
        this->data()[i] = *p;
        i++;
    }
}

/**
 * access the matrix element by index with bound check.
 *
 * @ param i the row index of the element
 * @ param j the column index of the element
 * @ return referent to the matrix element.
 */
const double & Matrix::at(size_t i, size_t j) const
{
    if (i >= row_ || j >= col_) {
        std::stringstream msg;
        msg << "Index is out of range at position (" << i << ", " << j << ".";
        throw exception::IndexRangeError(msg.str());
    }
    return (*this)(i, j);
}

/**
 * Initialize matrix object with the help of `<<` and `,` operator in an easy way.
 * That is, `A << 1, 2, 3, 4;`. The matrix dimension of A has to be specified
 * in advance and the number of elements has to be equal to the matrix size,
 * otherwise, it will abort with error.
 */
MatrixCommaInitializer Matrix::operator <<(double a)
{
    return MatrixCommaInitializer(*this, a);
}

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
            printf(" %15.8e,", this->data()[k]);
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
            printf(" %15.8e,", this->data()[ij]);
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

/**
 * check if two matrix has the same dimension.
 *
 * @ param[in] other: the other matrix.
 * @ return bool.
 */
bool Matrix::is_same_dimension_to(const Matrix &other) const
{
    return ((this->row() == other.row()) && (this->col() == other.col()));
}


/**
 * calculate matrix trace.
 */
double Matrix::trace() const
{
    if (!this->is_square()) {
        throw exception::DimensionError("Cannot get trace of a matrix that is not squared.");
    }
    double rst = 0.0;
    for (size_t i = 0; i < row_; i++) {
        rst += (*this)(i, i);
    }
    return rst;
}

/**
 * Make the matrix to be symmetric.
 *
 * @ param[in] uplo "U" using upper triangular part; "L" using lower triangular part.
 * @ return *this the symmetrized matrix itself.
 */
Matrix & Matrix::to_symmetric(const string & uplo)
{
    if (row_ != col_) {
        throw exception::DimensionError("Cannot symmetrize a matrix that is not squared.");
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
    blas::dscal_(&size, &alpha, this->data(), blas::ione);
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
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < size_; i++) {
        this->data()[i] = a;
    }
}

/**
 * set matrix to be identity.
 */
Matrix & Matrix::set_identity()
{
    if (!this->is_square()) {
        throw exception::DimensionError("Cannot make a non-square matrix to be identity.");
    }
    this->fill_all(0.0);
    #ifdef USE_OPENMP
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < row_; i++) {
        (*this)(i, i) = 1.0;
    }
}

}   // namespace matrix
