#include "matrix_io.h"
#include "matrix.h"
#include <unistd.h>
#include <assert.h>

namespace matrix {
/**
 * write a number of matrix into binary file in order.
 */
void write_matrices_to_binary(vector<std::shared_ptr<const Matrix>> &Mat, const char *fname)
{
    if (fname == NULL) {
        printf("Error in mtx_write_matrices(): file name is not given!\n");
        std::exit(EXIT_FAILURE);
    } else if(access(fname, F_OK) == 0) {
        if(access(fname, W_OK) != 0 ) {
            printf("Error in mtx_write_matrices(): file \"%s\" "
                   "has no writing access.\n", fname);
            std::exit(EXIT_FAILURE);
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

/**
 * read a number of matrices from a binary file.
 */
void read_matrices_from_binary(vector<std::shared_ptr<Matrix>> &Mat, const char *fname)
{
    assert(fname);
    if(access(fname, R_OK) != 0 ) {
        printf("Error in mtx_read_matrices(): file \"%s\" either does not exist"
               " or has no reading access.\n", fname);
        std::exit(EXIT_FAILURE);
    }

    FILE *f = fopen(fname, "rb");
    for (size_t i = 0; i < Mat.size(); i++) {
        size_t read_row = 0;
        size_t read_col = 0;
        fread(&read_row, sizeof(read_row), 1, f);
        fread(&read_col, sizeof(read_col), 1, f);
        if (read_row != Mat[i]->row() || read_col != Mat[i]->col()) {
            printf("Error in read %zu-th matrix from binary file \"%s\":"
                   "matrix dimensions do not match!", i, fname);
            std::exit(EXIT_FAILURE);
        }
        fread(Mat[i]->data(), sizeof (double), read_row * read_col, f);
    }
    fclose(f);
}

}
