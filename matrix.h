#pragma once
#include <stdlib.h>


typedef struct Matrix {
    float* data;
    size_t rows;
    size_t cols;
} Matrix;

Matrix* mCreate(size_t rows, size_t cols);
void mFree(Matrix* matrix);
Matrix* mCopy(Matrix* matrix);

float* mRow(Matrix* mat, size_t index);

void mPrint(Matrix* mat);

void mTranspose(Matrix* mat);
Matrix* mCreateTranspose(Matrix* mat);
void mScale(Matrix* mat, float scalar);
Matrix* mAdd(Matrix* a, Matrix* b);
Matrix* mMultiply(Matrix* a, Matrix* b);