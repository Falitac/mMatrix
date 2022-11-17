#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

static Matrix* mAlloc(size_t rows, size_t cols) {
    if(rows == 0 || cols == 0) {
        return NULL;
    }
    Matrix* result = malloc(sizeof(Matrix));
    if(!result) {
        return NULL;
    }
    result->data = malloc(sizeof(float) * rows * cols);
    if(!result->data) {
        free(result);
        return NULL;
    }
    result->rows = rows;
    result->cols = cols;

    return result;
}

Matrix* mCreate(size_t rows, size_t cols) {
    Matrix* result = mAlloc(rows, cols);
    return result;
}

void mFree(Matrix* matrix) {
    if(matrix == NULL) {
        return;
    }
    free(matrix->data);
    free(matrix);
    matrix = NULL;
}


Matrix* mCopy(Matrix* matrix) {
    if(!matrix) return NULL;
    Matrix* result = mAlloc(matrix->rows, matrix->cols);
    memcpy(result->data, matrix->data, sizeof(float) * matrix->rows * matrix->cols);
    return result;
}

void mFill(Matrix* mat, float val) {
    if(!mat) return;
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            mRow(mat, i)[j] = val;
        }
    }
}

void mGenerate(Matrix* mat, float (*generator)(Matrix* mat, int i, int j)) {
    if(!mat) return;
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            mRow(mat, i)[j] = generator(mat, i, j);
        }
    }
}

Matrix* mLoadFromFile(const char* filename) {
    printf("Loading from: %s\n", filename);
    FILE* file = fopen(filename, "r");
    if(!file) {
        return NULL;
    }

    size_t rows, cols;
    int scannedNumber = fscanf(file, "%zu %zu", &rows, &cols);
    if(scannedNumber != 2) {
        return NULL;
    }

    Matrix* result = mAlloc(rows, cols);
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            scannedNumber = fscanf(file, "%f", &mRow(result, i)[j]);
            if(scannedNumber != 1) {
                mFree(result);
                return NULL;
            }
        }
    }
    fclose(file);

    return result;
}

float* mRow(Matrix* mat, size_t index) {
    return &mat->data[index * mat->cols];
}

void mPrint(Matrix* mat) {
    if(!mat) return;
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            printf("%6.2f ", mRow(mat, i)[j]);
        }
        puts("");
    }
}

void mTranspose(Matrix* mat) {
    if(!mat) return;
    size_t tmp = mat->rows;
    mat->rows = mat->cols;
    mat->cols = tmp;
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
                mRow(result, i)[j] += mRow(a, i)[k] * mRow(b, k)[j];
            }
        }
    }

    return result;
}

/* TODO: This is temporary solution, looking for better one*/
void mMulFirst(Matrix* a, Matrix* b) {
    Matrix* result = mMul(a, b);
    mFree(a);
    *a = *result;
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

static int mFindMaxAbs(Matrix* mat, int index) {
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
    mSwapRows(mat, index, mFindMaxAbs(mat, index));
}

void mReduceRows(Matrix* mat) {
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

float mDiagonalProduct(Matrix* mat) {
    size_t lower = mat->rows < mat->cols ? mat->rows : mat->cols;
    float result = 1.;
    for(int i = 0; i < lower; i++) {
        result *= mRow(mat, i)[i];
    }
    return result;
}

float mDot(Matrix* vec1, Matrix* vec2) {
    assert(vec1->rows == vec2->rows && vec1->cols == vec2->cols);
    float accum = 0.f;
    for(int i = 0; i < vec1->rows; i++) {
        accum += mRow(vec1, 0)[i] * mRow(vec2, 0)[i];
    }
    return accum;
}

void mNormalize(Matrix* vec) {
    assert(vec->cols == 1);
    float accum = sqrt(mDot(vec, vec));
    mScale(vec, 1 / accum);
}