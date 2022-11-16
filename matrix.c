#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

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

void mGenerate(Matrix* mat, float (*generator)(Matrix* mat, int i, int j)) {
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            mRow(mat, i)[j] = generator(mat, i, j);
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

Matrix* mCreateMinor(Matrix* mat, int index1, int index2) {
    Matrix* result = mAlloc(mat->rows - 1, mat->cols - 1);
    float* it = result->data;

    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            if(i == index1 || j == index2) continue;
            *it++ = mRow(mat, i)[j];
        }
    }
    return result;
}

float mDet(Matrix* mat) {
    assert(mat->rows == mat->cols);
    if(mat->rows == 1) return mat->data[0];

    float result = 0.f;
    for(int i = 0; i < mat->cols; i++) {
        float sign = (i&1) ? -1.f : 1.f;
        
        Matrix* minorMat = mCreateMinor(mat, 0, i);
        result += sign * mRow(mat, 0)[i] * mDet(minorMat);
        mFree(minorMat);
    }
    return result;
}

void mSwapRows(Matrix* mat, int index1, int index2) {
    if(index1 == index2) return;
    for(int i = 0; i < mat->cols; i++) {
        float tmp = mRow(mat, index1)[i];
        mRow(mat, index1)[i] = mRow(mat, index2)[i];
        mRow(mat, index2)[i] = tmp;
    }
}

static int mFindMax(Matrix* mat, int index) {
    float maxAbs = fabs(mRow(mat, index)[index]);
    int maxIndex = index;
    for(int i = index + 1; i < mat->rows; i++) {
        if(fabs(mRow(mat, i)[index]) > maxAbs) {
            maxIndex = i;
            maxAbs = fabs(mRow(mat, i)[index]);
        }
    }
    return maxIndex;
}

static void mPivotRows(Matrix* mat, int index) {
    mSwapRows(mat, index, mFindMax(mat, index));
}

void mRedRows(Matrix* mat) {
    for(int i = 0; i < mat->rows - 1; i++) {
        mPivotRows(mat, i);
        for(int j = i + 1; j < mat->cols; j++) {
            float coeff = mRow(mat, j)[i] / mRow(mat, i)[i];
            for(int k = i; k < mat->cols; k++) {
                mRow(mat, j)[k] -= coeff * mRow(mat, i)[k];
            }
        }
    }
}

float mFactorDiagonal(Matrix* mat) {
    size_t lower = mat->rows < mat->cols ? mat->rows : mat->cols;
    float result = 1.;
    for(int i = 0; i < lower; i++) {
        result *= mRow(mat, i)[i];
    }
    return result;
}