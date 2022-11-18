#pragma once
#include <stdlib.h>


typedef struct Matrix {
    float* data;
    size_t rows;
    size_t cols;
} Matrix;

Matrix* mCreate(size_t rows, size_t cols);
Matrix* mIdentity(size_t size);
void mFree(Matrix* matrix);
Matrix* mCopy(Matrix* matrix);
void mFill(Matrix* mat, float val);
void mGenerate(Matrix* mat, float (*generator)(Matrix* mat, int i, int j));
Matrix* mLoadFromFile(const char* filename);

float* mRow(Matrix* mat, size_t index);

void mPrint(Matrix* mat);

int mIsClose(Matrix* m1, Matrix* m2);

void mTranspose(Matrix** mat);
Matrix* mCreateTranspose(Matrix* mat);

void mScale(Matrix* mat, float scalar);
void mAdd(Matrix* a, Matrix* b);
Matrix* mMul(Matrix* a, Matrix* b);
void mMulFirst(Matrix** a, Matrix* b);

Matrix* mCreateMinor(Matrix* mat, int index1, int index2);
float mDet(Matrix* mat);

void mReduceRows(Matrix* mat);
float mDiagonalProduct(Matrix* mat);
float mTrace(Matrix* mat);

float mDot(Matrix* vec1, Matrix* vec2);
void mNormalize(Matrix* vec);

float mCalcMaxEigenSymm(Matrix* mat);
Matrix* mCalcEigVecSymm(Matrix* mat, float maxEig);


Matrix* mCalcEigVecSymmetric(Matrix* mat);
Matrix* mCalcEigensSymmetric3(Matrix* mat);

void mSVD3(Matrix* mat, Matrix** u, Matrix** s, Matrix** v);