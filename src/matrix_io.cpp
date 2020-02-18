/**
 * @file
 * @brief definition relates to matrix I/O.
 */
#include "matrix_io.h"

#include "exception.h"
#include "matrix.h"
#include <assert.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
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

/**
 * @note The txt file `fname` will always be overwritten if `Mat` is not empty.
 */
void write_matrices_to_txt(vector<std::shared_ptr<Matrix>> &Mat,
                           const string &fname, size_t num_per_line)
{
    if (Mat.size() == 0)
        return;

    FILE *f = fopen(fname.c_str(), "w");
    if (f == NULL)
        throw exception::MatrixIOException(
            fname, "Cannot open file to write matrices.");
    for (size_t i = 0; i < Mat.size(); ++i) {
        const Matrix &A = *(Mat[i]);
        if (i == 0)
            fprintf(f, "Dimension,%zu,%zu\n", A.row(), A.col());
        else
            fprintf(f, "\nDimension,%zu,%zu\n", A.row(), A.col());

        // loop over each elements.
        size_t n = 0;
        const size_t size = A.size();
        for (size_t ii = 0; ii < A.size(); ++ii) {
            const auto val = A.data()[ii];
            if (n != num_per_line && ii != size - 1) {
                fprintf(f, "%.16e,", val);
                ++n;
            } else if (n == num_per_line && ii != size - 1) {
                fprintf(f, "%.16e\n", val);
                n = 0;
            } else {
                fprintf(f, "%.16e", val);
            }
        }
    }
    fclose(f);
}

/**
 * @brief split string with given delimeter.
 * @param [in] str: input string to be splitted.
 * @param [in] delim: delimeter.
 * @param [in, out] str_split: On exit, it stores the splitted string.
 */
static void str_split(std::string &str, const char delim,
                      std::vector<std::string> &str_split)
{
    str_split.clear();
    std::stringstream str_stream(str + std::string(1, delim));
    std::string word;
    while (std::getline(str_stream, word, delim)) {
        str_split.push_back(word);
    }
}

std::vector<std::shared_ptr<Matrix>> read_matrices_from_txt(const string &fname)
{
    std::ifstream fin;
    fin.open(fname);
    if (!fin)
        throw exception::MatrixIOException(
            fname, "Cannot open file to read matrices.");

    bool is_first_line = true;
    std::string line;
    std::vector<std::string> line_split;
    std::vector<std::shared_ptr<Matrix>> rst;
    std::shared_ptr<Matrix> current_matrix;
    std::size_t n_row = 0;
    std::size_t n_col = 0;
    std::size_t count = 0;

    // loop over each line in the txt file
    while (std::getline(fin, line)) {
        str_split(line, ',', line_split);
        auto p_data = line_split.begin();
        if (line_split[0] == "Dimension") {
            // now it starts to read a new matrix.
            // Before that, check if the last matrix is read successfully or
            // not.
            if (!is_first_line) {
                if (count != current_matrix->size()) {
                    throw exception::MatrixIOException(
                        fname, "Error in read matrix, unmatched element size.");
                }
            }
            // read dimension information for the new matrix.
            if (line_split.size() < 3) {
                throw exception::MatrixIOException(
                    fname, "Cannot read matrix dimension.");
            }
            try {
                n_row = std::stoi(line_split[1]);
                n_col = std::stoi(line_split[2]);
            } catch (...) {
                throw exception::MatrixIOException(
                    fname, "Failed to read matrix dimension.");
            }
            // create a new matrix and push it back to matrix vector.
            current_matrix = std::make_shared<Matrix>(n_row, n_col);
            rst.push_back(current_matrix);
            p_data += 3;
            count = 0; // set matrix element count as zero.
        }
        // put matrix element in current line into the matrix object.
        for (; p_data != line_split.end(); ++p_data) {
            try {
                current_matrix->data()[count] = std::stof(*p_data);
                ++count;
            } catch (...) {
                throw exception::MatrixIOException(
                    fname, "Failed to read matrix element.");
            }
        }
        // at this point, the first line has been parsed.
        is_first_line = false;
    }
    fin.close();
    return rst;
}

} // namespace matrix
