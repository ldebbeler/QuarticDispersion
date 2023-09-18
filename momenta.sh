#!/bin/bash

g++ -std=c++17 -Wall -g -pedantic -I ./include  \
    src/interpolate1d.cpp                       \
    src/hcubature.c                             \
    src/scaling.cpp                             \
    src/grid.cpp                                \
    src/bubble.cpp                              \
    src/selfenergy.cpp                          \
    src/testBubble.cpp                          \
    src/writeReal.cpp                           \
    src/momentumDep.cpp                         \
    -o bin/MomentumDep.elf  -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl

echo "Compilation finished!"

#./bin/MomentumDep.elf
