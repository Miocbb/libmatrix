/**
 * file: io.hpp
 *
 * I/O related.
 */

#ifndef _MATRIX_IO_H_
#define _MATRIX_IO_H_

#include <iostream>
#include <cstdlib>

namespace matrix {

/**
 * print out error message and exit program.
 */
inline void sig_err(char *msg)
{
    std::cout << msg << "\n";
    std::exit(EXIT_FAILURE);
}
}

#endif
