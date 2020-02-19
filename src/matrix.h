/**
 * @file
 * @brief Decleration matrix class.
 */

#ifndef _MATRIX_SRC_MATRIX_H_
#define _MATRIX_SRC_MATRIX_H_

#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

/**
 * @brief top-level matrix namespace.
 */
namespace matrix {

using std::size_t;
using std::string;
using std::vector;

/**
 * @brief Helper class used to do comma initialization for matrix object like
 * Eigen library.
 */
class MatrixCommaInitializer;

/**
 * @brief Matrix class declaration.
 */
class Matrix {
  private:
    size_t row_;
    size_t col_;
    size_t size_;
    /* The data of matrix is either stored in vector `data_` or outside of the
     * object that is pointed by the double pointer `data_ptr_`. */
    vector<double> data_vec_;
    double *data_ptr_; /* This pointer will always point to the head of the
                        * matrix data.
                        * No memory allocation is assigned for this pointer. */

  public:
    /**
     * @brief copy types.
     */
    enum CopyType {
        kShallowCopy, /**< shallow copy: only gain data access and avoid copying
                         data elements */
        kDeepCopy,    /**< deep copy: copy all the data. */
    };

    /**
     * @brief Construct a matrix with given size.
     *
     * @details The memory for a matrix with input size will
     * be allocated and initialize all matrix elements
     * to be zero. Matrix data is stored in inside of this object.
     *
     * @param [in] row: Number of rows of the matrix
     * @param [in] col: Number of columns of the matrix
     */
    Matrix(size_t row, size_t col)
        : row_(row), col_(col), data_vec_(row * col, 0.0),
          size_(row * col), data_ptr_{data_vec_.data()}
    {
    }

    /**
     * @brief Construct a matrix from a std::vector<double>.
     *
     * @param [in] row: Number of rows of the matrix.
     * @param [in] col: Number of columns of the matrix.
     * @param [in] inp_data: A std::vector<double> that stores the matrix data.
     * @param [in] copy_type: Copy type to initialize the matrix from \p
     * inp_data. If \p copy_type equals to matrix::Matrix::kDeepCopy, all the
     * data will be copied into matrix. If the \p copy type equals to
     * matrix::Matrix::kShallowCopy, the pointer to the data memory is assigned
     * to the matrix and the matrix gains the access to the data. Default is
     * deep copy.
     *
     * @note The data size of the vector will be check to match the matrix size.
     * @note This constructor is not const friendly. If you want to do a shallow
     * copy of a const vector to the matrix, you can do declare the created
     * Matrix to be const and typecast the const vector to be non-const vector,
     * that is,
     * @code
     * const vector<double> data;
     * const Matrix A(row, col, const_cast <vector<double>>(data),
     *                kShallowCopy);
     * @endcode
     * This can preserve the created matrix not change the original vector
     * data.
     */
    Matrix(size_t row, size_t col, vector<double> &inp_data,
           CopyType copy_type = kDeepCopy);

    /**
     * @brief Construct a matrix from a double array pointer.
     *
     * @param [in] inp_data_ptr: a pointer points to double array that stores
     * the matrix data.
     * @param [in] copy_type: If \p copy_type equals matrix::Matrix::kDeepCopy,
     * all the data will be copied into matrix. If the \p copy_type equals to
     * matrix::Matrix::kShallowCopy, the pointer \p inp_data_ptr is stored in
     * the matrix to gain the access to matrix data. Default is deep copy.
     *
     * @note The data size of the array will NOT be check to match the matrix
     * size. Use this constructor carefully.
     * @note This constructor is not const friendly. If you want to do a shallow
     * copy of a const double array to the matrix, you can do declare the
     * created Matrix to be const and typecast the const array to be non-const
     * array, that is,
     * ```
     * const double * data;
     * const Matrix A(row, col, const_cast <double *>(data), kShallowCopy);
     * ```
     * This can preserve the created matrix not change the original array data.
     */
    Matrix(size_t row, size_t col, double *inp_data_ptr,
           CopyType copy_type = kDeepCopy);

    /**
     * @brief Default constructor.
     * @details Create an empty matrix object with dimension [0, 0].
     */
    Matrix() : row_(0), col_(0), size_(0), data_vec_(0), data_ptr_{nullptr} {}

    /**
     * @brief Copy constructor: copy from a matrix.
     *
     * @param [in] other: the other matrix to be copied.
     */
    Matrix(const Matrix &other);

    /**
     * @brief Copy assignment operator: assign from a matrix.
     *
     * @param [in] other: the other matrix used for assignment.
     */
    Matrix &operator=(const Matrix &other);

    /**
     * @brief Copy assignment operator: enable an easy way to do matrix
     * element assignment from std::initializer_list.
     *
     * @ param[in] init_list: the data contained in the initializer_list.
     * @ return: the const reference to the matrix itself.
     *
     * @note The matrix dimension has to be declared in advance. The input list
     * size will be checked for the initialization. Below is an example,
     * @code
     * Matrix A(2, 2); // create a matrix with dimension [2, 2] first.
     * A = {1, 2,      // initialization with initializer_list. A(0, 1) = 2,
     *      3, 4};     // and A(1, 0) = 3.
     * A = {1, 2, 3};  // Error: will throw an error for unmatched size.
     * @endcode
     */
    const Matrix &operator=(std::initializer_list<double> init_list);

    /**
     * @brief operator << overloading for comma initialization like Eigen3
     * library.
     *
     * @param [in] a: an element.
     *
     * @note Initialize matrix object with the help of operator `<<` and `,` in
     * an easy way. That is,
     * @code
     * Matrix A(2, 2); // create a matrix with dimension [2, 2] first.
     * A << 1, 2,      // initialization with comma initialization. A(0, 1) = 2,
     *      3, 4;      // and A(1, 0) = 3.
     * A << 1, 2, 3;   // Error: will throw an error for unmatched size.
     * @endcode
     * The matrix dimension of A has to be specified in advance and the
     * number of elements has to be equal to the matrix size, otherwise, it
     * will throw an error.
     */
    MatrixCommaInitializer operator<<(double a);

    /**
     * @brief Access/modify the matrix element by index without
     * bound check.
     *
     * @param [in] i: The row index of the element
     * @param [in] j: The column index of the element
     * @return double&: matrix element [\p i, \p j].
     */
    double &operator()(size_t i, size_t j)
    {
        return const_cast<double &>(static_cast<const Matrix &>(*this)(i, j));
    };

    /**
     * @brief Access the matrix element by index without bound check.
     *
     * @param [in] i: The row index of the element
     * @param [in] j: The column index of the element
     * @return const double&: matrix element [\p i, \p j].
     */
    const double &operator()(size_t i, size_t j) const
    {
        return data_ptr_[i * col_ + j];
    };

    /**
     * @brief Access/modify the matrix element by index with bound check.
     *
     * @param [in] i: The row index of the element
     * @param [in] j: The column index of the element
     * @return double&: matrix element [\p i, \p j].
     */
    double &at(size_t i, size_t j)
    {
        return const_cast<double &>(
            static_cast<const Matrix &>(*this).at(i, j));
    }

    /**
     * @brief Access the matrix element by index with bound check.
     *
     * @param [in] i: The row index of the element
     * @param [in] j: The column index of the element
     * @return const double&: matrix element [\p i, \p j].
     */
    const double &at(size_t i, size_t j) const;

    /**
     * @brief Get the pointer that points to the begining of matrix data.
     * @return double *
     */
    double *data() { return data_ptr_; }

    /**
     * @brief Get the const pointer that points to the begining of matrix data.
     * @return const double *
     */
    const double *data() const { return data_ptr_; }

    /**
     * @brief Get matrix row numbers.
     * @return size_t
     */
    const size_t &row() const { return row_; }

    /**
     * @brief Get matrix column numbers.
     * @return size_t
     */
    const size_t &col() const { return col_; }

    /**
     * @brief Get matrix size, that is, the total number of matrix elements.
     * @return size_t
     */
    const size_t &size() const { return size_; }

    /**
     * @brief Check if the matrix data is managed by the object.
     * @details If true, the matrix data will be destroyed along with the
     * object. otherwise, it will not.
     * @return bool
     * @note if the matrix is an empty matrix, it will return false.
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
     * @brief Check if the matrix is square or not.
     * @return bool
     */
    bool is_square() const { return (row_ == col_); }

    /**
     * @brief Check if the current matrix is symmetric or not based on
     * the input threshold. Default threshold is 1e-10.
     *
     * @param [in] threshold: the threshold of testing equality between two
     * float numbers.
     * @return bool.
     */
    bool is_symmetric(double threshold = 1e-10) const;

    /**
     * @brief Check if the current matrix is a diagonal matrix or not based on
     * input threshold. Default threshold is 1e-10.
     *
     * @param [in] threshold: the threshold of testing equality between two
     * float numbers.
     * @return bool.
     */
    bool is_diagonal(double threshold = 1e-10) const;

    /**
     * @brief Check if the current matrix is identity or not based on
     * input threshold. Default threshold is 1e-10.
     *
     * @param [in] threshold: the threshold of testing equality between two
     * float numbers.
     * @return bool.
     */
    bool is_identity(double threshold = 1e-10) const;

    /**
     * @brief Check if the current matrix is a zero matrix or not based on
     * input threshold. Default threshold is 1e-10.
     *
     * @param [in] threshold: the threshold of testing equality between two
     * float numbers.
     * @return bool.
     */
    bool is_zeros(double threshold = 1e-10) const;

    /**
     * @brief Check if two matrix is equal by scanning everything based on
     * a threshold. By default, the threshold is 1e-10.
     *
     * @param [in] other: the other matrix to be compared with.
     * @return bool
     */
    bool is_equal_to(const Matrix &other, double threshold = 1e-10) const;

    /**
     * @brief Check if two matrix has the same dimension.
     *
     * @param[in] other: the other matrix to be compared with.
     * @return bool.
     */
    bool is_same_dimension_to(const Matrix &other) const;

    /**
     * @brief print out the full matrix.
     * @param [in] elements_per_line: number of elements being printed per line.
     *  Default is 5.
     */
    void show_full(size_t elements_per_line = 5) const;

    /**
     * @brief Print out the lower triangular matrix (including diagonal
     * elements).
     * @param [in] elements_per_line: number of elements being printed per line.
     *  Default is 5.
     */
    void show_lower(size_t elements_per_line = 5) const;

    /**
     * @brief Calculate matrix trace.
     * @return double: the matrix trace.
     */
    double trace() const;

    /**
     * @brief Resize the matrix into given dimension.
     * @details Matrix data value is perserved. The order of data stored in
     * memory is preserved as well, so you will need to re-interprate the new
     * matrix based on the new size.
     *
     * When new matrix size is greater than the original matrix size,
     * zeros will be appended to the end.
     *
     * @param [in] row: new number of rows.
     * @param [in] col: new number of columns.
     */
    void resize(size_t row, size_t col);

    /**
     * @brief Make the matrix to be symmetric.
     *
     * @param[in] uplo: when \p uplo equals to "U", the upper triangular part is
     * used. when \p uplo equals to "L", the lower triangular part is used.
     */
    void to_symmetric(const string &uplo);

    /**
     * @brief Make the matrix to be random with elements uniformly distributed
     * in range [a, b).
     *
     * @param [in] a: left range bound.
     * @param [in] b: right range bound.
     *
     * @note The random number generator is initialized with a NON-FIXED seed.
     * So the randomness behavior is not repeatable at running time.
     */
    void randomize(double a, double b);

    /**
     * @brief Make the matrix to be random with elements uniformly distributed
     * in range [a, b).
     *
     * @param [in] a: left range bound.
     * @param [in] b: right range bound.
     *
     * @note The random number generator is initialized with a FIXED seed. So
     * the randomness behavior is repeatable at running time.
     */
    void randomize_seed_fixed(double a, double b);

    /**
     * @brief Scales current matrix by a constant.
     * @details A = alpha * A.
     * @param [in] alpha: the scalar coefficient.
     */
    void scale(const double alpha);

    /**
     * @brief Fill all elements with input number.
     * @details For all matrix element A(i, i), make A(i, i) = \p a
     * @param [in] a: the number to be filled.
     */
    void fill_all(double a);

    /**
     * @brief Set the matrix to be identity.
     */
    void set_identity();

    /**
     * @brief Set the matrix to be its transpose.
     */
    void transpose();
};

} // namespace matrix

#endif // _MATRIX_SRC_MATRIX_H_
