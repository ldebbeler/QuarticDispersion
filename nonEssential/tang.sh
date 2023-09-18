#!/bin/bash

g++ -std=c++17 -Wall -g -pedantic -I ./include  \
    src/interpolate1d.cpp                       \
    src/hcubature.c                             \
    src/scaling.cpp                             \
    src/grid.cpp                                \
    src/bubble.cpp                              \
    src/testBubble.cpp                          \
    src/writeTang.cpp                           \
    src/tangFlow.cpp                            \
    src/mainFlow.cpp                            \
    -o bin/mainFlow.elf  -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl

echo "Compilation finished!"

#./bin/mainFlow.elf
