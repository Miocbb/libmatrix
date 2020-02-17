/**
 * @file
 * @brief definition relates to exceptions used in matrix library.
 */
#ifndef _MATRIX_SRC_EXCEPTION_H_
#include "exception.h"

#include <sstream> // std::stringstream

namespace matrix {
namespace exception {

void MatrixException::make_message(const char *msg)
{
    msg_ << std::endl;
    msg_ << "Fatal error: " << msg << std::endl;
}

} // namespace exception
} // namespace matrix

#endif // _MATRIX_SRC_EXCEPTION_H_
