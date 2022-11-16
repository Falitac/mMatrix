#pragma once
#include <stdlib.h>

typedef struct Matrix {
    float* data;
    size_t rows;
    size_t cols;
} Matrix;

Matrix* genMatrix(size_t rows, size_t cols) {
    if(rows == 0 || cols == 0) {
        return NULL;
    }
}

void freeMatrix(Matrix* matrix) {
    free(matrix->data);
    free(matrix);
}