#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

static Matrix* mAlloc(size_t rows, size_t cols) {
    Matrix* result = malloc(sizeof(Matrix));
    result->data = malloc(sizeof(float) * rows * cols);
    result->rows = rows;
    result->cols = cols;
}

Matrix* mCreate(size_t rows, size_t cols) {
    if(rows == 0 || cols == 0) {
        return NULL;
    }
    
    Matrix* result = mAlloc(rows, cols);
    return result;
}

void mFree(Matrix* matrix) {
    if(matrix == NULL) {
        return;
    }
    free(matrix->data);
    free(matrix);
}


Matrix* mCopy(Matrix* matrix) {
    Matrix* result = mAlloc(matrix->rows, matrix->cols);
    memcpy(result->data, matrix->data, sizeof(float) * matrix->rows * matrix->cols);
    return result;
}

void mFill(Matrix* mat, float val) {
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            mRow(mat, i)[j] = val;
        }
    }
}

float* mRow(Matrix* mat, size_t index) {
    return &mat->data[index * mat->cols];
}

void mPrint(Matrix* mat) {
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            printf("%5.3f ", mRow(mat, i)[j]);
        }
        puts("");
    }
}

void mTranspose(Matrix* mat) {
    for(int i = 0; i < mat->rows; i++) {
        for(int j = i + 1; j < mat->cols; j++) {
            float tmp = mRow(mat, i)[j];
            mRow(mat, i)[j] = mRow(mat, j)[i];
            mRow(mat, j)[i] = tmp;
        }
    }
}

Matrix* mCreateTranspose(Matrix* mat) {
    Matrix* result = mCopy(mat);
    mTranspose(result);
    return result;
}

void mScale(Matrix* mat, float scalar) {
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            mRow(mat, j)[i] *= scalar;
        }
    }
}

void mAdd(Matrix* a, Matrix* b) {
    assert(a->cols == b->cols && a->rows == b->cols);
    for(int i = 0; i < a->rows; i++) {
        for(int j = 0; j < a->cols; j++) {
            mRow(a, i)[j] += mRow(b, i)[j];
        }
    }
}

Matrix* mMul(Matrix* a, Matrix* b) {
    assert(a->cols == b->rows);
    Matrix* result = mAlloc(a->rows, b->cols);
    mFill(result, 0.f);
    for(int i = 0; i < result->rows; i++) {
        for(int j = 0; j < result->cols; j++) {
            for(int k = 0; k < a->rows; k++) {
                mRow(result, i)[j] = mRow(a, i)[k] * mRow(b, k)[j];
            }
        }
    }

    return result;
}