#ifndef CONST_H_
 #define CONST_H_

 #define NUM_OF_SYMBOLS 3  // Impacts vertex.* files, stands for s0-s2 symbols
 #define ORDER_CONSTANT 2  // Controls the maximum size of the primary heap
 #define EMPTY_LABEL "__PHI__"  // Impacts vertex.cpp, computeChiSqValue

 #define TOPK 10  // Defines the highest valued chi squared vertex candidates considered in the greedy approach

 #define LAPLACIAN_BIAS 0.0001 // Corrective laplacian bias for chi-square, useful when the expected symbol occurrence value is very small

#endif
