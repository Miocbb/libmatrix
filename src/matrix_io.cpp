/**
 * @file
 * @brief definition relates to matrix I/O.
 */
#include "matrix_io.h"

#include "exception.h"
#include "matrix.h"
#include <assert.h>
#include <sstream>
#include <unistd.h>

namespace matrix {
void write_matrices_to_binary(vector<std::shared_ptr<const Matrix>> &Mat,
                              const char *fname)
{
    if (fname == NULL) {
        throw exception::MatrixIOException(fname, "No file name.");
    } else if (access(fname, F_OK) == 0) {
        if (access(fname, W_OK) != 0) {
            throw exception::MatrixIOException(fname, "No wrting access.");
        }
    }

    FILE *f = fopen(fname, "wb");
    for (size_t i = 0; i < Mat.size(); i++) {
        size_t row = Mat[i]->row();
        size_t col = Mat[i]->col();
        fwrite(&row, sizeof(row), 1, f);
        fwrite(&col, sizeof(col), 1, f);
        fwrite(Mat[i]->data(), sizeof(double), row * col, f);
    }
    fclose(f);
}

void read_matrices_from_binary(vector<std::shared_ptr<Matrix>> &Mat,
                               const char *fname)
{
    if (fname == NULL) {
        throw exception::MatrixIOException(fname, "No file name.");
    } else if (access(fname, R_OK) != 0) {
        throw exception::MatrixIOException(fname, "No readding access.");
    }

    FILE *f = fopen(fname, "rb");
    for (size_t i = 0; i < Mat.size(); i++) {
        size_t read_row = 0;
        size_t read_col = 0;
        fread(&read_row, sizeof(read_row), 1, f);
        fread(&read_col, sizeof(read_col), 1, f);
        if (read_row != Mat[i]->row() || read_col != Mat[i]->col()) {
            std::stringstream msg;
            msg << "Error in read " << i << "-th matrix from binary file"
                << fname << ". Dimension is not matched.";
            throw exception::MatrixIOException(fname, msg.str());
        }
        fread(Mat[i]->data(), sizeof(double), read_row * read_col, f);
    }
    fclose(f);
}

std::vector<std::shared_ptr<Matrix>>
read_matrices_from_binary(const char *fname)
{
    if (fname == NULL) {
        throw exception::MatrixIOException(fname, "No file name.");
    } else if (access(fname, R_OK) != 0) {
        throw exception::MatrixIOException(fname, "No readding access.");
    }

    std::vector<std::shared_ptr<Matrix>> rst;
    FILE *f = fopen(fname, "rb");
    while (!feof(f)) {
        size_t read_row = 0;
        size_t read_col = 0;
        fread(&read_row, sizeof(read_row), 1, f);
        fread(&read_col, sizeof(read_col), 1, f);

        auto mat = std::make_shared<Matrix>(read_row, read_col);
        fread(mat->data(), sizeof(double), read_row * read_col, f);
        rst.push_back(mat);
    }
    fclose(f);

    return rst;
}

std::vector<std::shared_ptr<Matrix>> read_matrices_from_binary(string &fname)
{
    return read_matrices_from_binary(fname.c_str());
}

} // namespace matrix
