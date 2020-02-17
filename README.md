# libmatrix

This is a library that wraps the commonly used blas/lapack subroutines for
matrix operations. Use this library to make life easier.

How to compile:
c++11 standard is required.
```
mkdir build
cd build
cmake ../
make
```

To generate documentation: need `doxygen` tool.
```
cd doc
doxygen Doxyfile
```
Then open the `doc/doxygen_html/index.html` with your preferred browser to check it out.
