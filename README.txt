This directory contains the code used to compute the results presented in Phys Rev B 107, 165152 (2023)
Contains:
- arXiv preprint of the article, contains the same results
- include directory: header files and files to facilitate use of external libraries
- src directory: source files to corresponding headers. main function in 'mainCalculation.cpp'
- data directory: jupyter notebooks to visualize results


Prerequisites:
- C++17 standard
- hdf5 library: H5CPP and H5PY (C++ and Python)
- GSL library

Usage:
*) create a directory "mkdir bin". The executable file is going to be stored in this directory.
*) Make sure "exec.sh" is executable and run ./exec.sh
*) This involves compiling, linking and execution
*) .h5 file is generated. Name can be changed in "include/constants.h"
*) file can be loaded into jupyter notebook and read with the provided class
*) The class in the jupyter notebook is able to provide all graphical results
*) calculation for static bare susceptibility based on tight binding model is provided in another jupyter notebook

