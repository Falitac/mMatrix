#pragma once
#include <stdlib.h>


typedef struct Matrix {
    double* data;
    size_t rows;
    size_t cols;
} Matrix;

void mSetPrintNumberFormat(const char* str);
Matrix* mCreate(size_t rows, size_t cols);
Matrix* mIdentity(size_t size);
void mFree(Matrix* matrix);
Matrix* mCopy(Matrix* matrix);
void mFill(Matrix* mat, double val);
void mGenerate(Matrix* mat, double (*generator)(Matrix* mat, int i, int j));
Matrix* mLoadFromFile(const char* filename);

double* mRow(Matrix* mat, size_t index);

void mPrint(Matrix* mat);

int mIsClose(Matrix* m1, Matrix* m2);

void mTranspose(Matrix** mat);
Matrix* mCreateTranspose(Matrix* mat);

void mScale(Matrix* mat, double scalar);
void mAdd(Matrix* a, Matrix* b);
Matrix* mMul(Matrix* a, Matrix* b);
void mMulFirst(Matrix** a, Matrix* b);

Matrix* mCreateMinor(Matrix* mat, int index1, int index2);
double mDet(Matrix* mat);

void mReduceRows(Matrix* mat);
double mDiagonalProduct(Matrix* mat);
double mTrace(Matrix* mat);

double mDot(Matrix* vec1, Matrix* vec2);
void mNormalize(Matrix* vec);

double mCalcMaxEigenSymm(Matrix* mat);
Matrix* mCalcEigVecSymm(Matrix* mat, double maxEig);
Matrix* mCalcEquivalentMinorMatrix(Matrix* mat, Matrix* x);

void mSVD(Matrix* mat, Matrix** u, Matrix** s, Matrix** v);

// this might be specialized
Matrix* mCalcEigVecSymmIterative(Matrix* mat);
Matrix* mCalcEigensSymmetric3Cardano(Matrix* mat);

void mSVD3(Matrix* mat, Matrix** u, Matrix** s, Matrix** v);