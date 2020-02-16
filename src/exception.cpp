#include "exception.h"

#include <sstream>  // std::stringstream

namespace matrix {
namespace exception {

void MatrixException::make_message(const char* msg)
{
    msg_ << std::endl;
    msg_ << "Fatal error: " << msg << std::endl;
}

}   // namespace matrix::exception
}   // namespace matrix
