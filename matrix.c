#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

static char G_M_FORMAT_BUFFER[64] = "%6.3f ";

void mSetPrintNumberFormat(const char* str) {
    size_t len = strlen(str);
    if(len >= 64) {
        return;
    }
    strncpy(G_M_FORMAT_BUFFER, str, len);
}

static void swap(double *xp, double *yp) {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}
 
static Matrix* mAlloc(size_t rows, size_t cols) {
    if(rows == 0 || cols == 0) {
        return NULL;
    }
    Matrix* result = malloc(sizeof(Matrix));
    if(!result) {
        return NULL;
    }
    result->data = malloc(sizeof(double) * rows * cols);
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
    mFill(result, 0.);
    return result;
}

Matrix* mIdentity(size_t size) {
    Matrix* result = mCreate(size, size);
    for(size_t i = 0; i < size; i++) {
        mRow(result, i)[i] = 1.;
    }
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
    memcpy(result->data, matrix->data, sizeof(double) * matrix->rows * matrix->cols);
    return result;
}

void mFill(Matrix* mat, double val) {
    if(!mat) return;
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            mRow(mat, i)[j] = val;
        }
    }
}

void mGenerate(Matrix* mat, double (*generator)(Matrix* mat, int i, int j)) {
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

double* mRow(Matrix* mat, size_t index) {
    return &mat->data[index * mat->cols];
}

void mPrint(Matrix* mat) {
    if(!mat) return;
    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            printf(G_M_FORMAT_BUFFER, mRow(mat, i)[j]);
        }
        puts("");
    }
}

int mIsClose(Matrix* m1, Matrix* m2) {
    assert(m1->cols == m2->cols && m1->rows == m2->rows);
    const double EPS = 1e-7;
    for(int i = 0; i < m1->cols * m1->cols; i++) {
        if(fabs(m1->data[i] - m2->data[i]) > EPS) {
            return 0;
        }
    }
    return 1;
}

void mTranspose(Matrix** mat) {
    if(!mat) return;
    Matrix* result = mCreateTranspose(*mat);
    mFree(*mat);
    *mat = result;
}

Matrix* mCreateTranspose(Matrix* mat) {
    if(!mat) return NULL;
    Matrix* result = mAlloc(mat->cols, mat->rows);

    for(int i = 0; i < result->rows; i++) {
        for(int j = 0; j < result->cols; j++) {
            mRow(result, i)[j] = mRow(mat, j)[i];
        }
    }
    return result;
}

void mScale(Matrix* mat, double scalar) {
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
            for(int k = 0; k < a->cols; k++) {
                mRow(result, i)[j] += mRow(a, i)[k] * mRow(b, k)[j];
            }
        }
    }

    return result;
}

/* TODO: This is temporary solution, looking for better one*/
void mMulFirst(Matrix** a, Matrix* b) {
    Matrix* result = mMul(*a, b);
    mFree(*a);
    *a = result;
}

Matrix* mCreateMinor(Matrix* mat, int index1, int index2) {
    Matrix* result = mAlloc(mat->rows - 1, mat->cols - 1);
    double* it = result->data;

    for(int i = 0; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            if(i == index1 || j == index2) continue;
            *it++ = mRow(mat, i)[j];
        }
    }
    return result;
}

double mDet(Matrix* mat) {
    assert(mat->rows == mat->cols);
    if(mat->rows == 1) return mat->data[0];

    double result = 0.f;
    for(int i = 0; i < mat->cols; i++) {
        double sign = (i&1) ? -1.f : 1.f;
        
        Matrix* minorMat = mCreateMinor(mat, 0, i);
        result += sign * mRow(mat, 0)[i] * mDet(minorMat);
        mFree(minorMat);
    }
    return result;
}

void mSwapRows(Matrix* mat, int index1, int index2) {
    if(index1 == index2) return;
    for(int i = 0; i < mat->cols; i++) {
        double tmp = mRow(mat, index1)[i];
        mRow(mat, index1)[i] = mRow(mat, index2)[i];
        mRow(mat, index2)[i] = tmp;
    }
}

static int mFindMaxAbs(Matrix* mat, int index) {
    double maxAbs = fabs(mRow(mat, index)[index]);
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
            double coeff = mRow(mat, j)[i] / mRow(mat, i)[i];
            for(int k = i; k < mat->cols; k++) {
                mRow(mat, j)[k] -= coeff * mRow(mat, i)[k];
            }
        }
    }
}

double mDiagonalProduct(Matrix* mat) {
    size_t lower = mat->rows < mat->cols ? mat->rows : mat->cols;
    double result = 1.;
    for(int i = 0; i < lower; i++) {
        result *= mRow(mat, i)[i];
    }
    return result;
}

double mTrace(Matrix* mat) {
    double result = 0.f;

    int lower = mat->rows < mat->cols ?  mat->rows : mat->cols;
    for(int i = 0; i < lower; i++) {
        result += mRow(mat, i)[i];
    }
    return result;
}

double mDot(Matrix* vec1, Matrix* vec2) {
    assert(vec1->rows == vec2->rows && vec1->cols == vec2->cols);
    double accum = 0.f;
    for(int i = 0; i < vec1->rows; i++) {
        accum += mRow(vec1, 0)[i] * mRow(vec2, 0)[i];
    }
    return accum;
}

void mNormalize(Matrix* vec) {
    assert(vec->cols == 1);
    double accum = sqrt(mDot(vec, vec));
    mScale(vec, 1 / accum);
}

double mCalcMaxEigenSymm(Matrix* mat) {
    Matrix* result = mAlloc(mat->rows, 1);
    mFill(result, 1.);

    Matrix* v = mMul(mat, result);
    double prevScalar = mRow(v, 0)[0] / mRow(result, 0)[0];
    for(int i = 0; i < 20; i++) {
        Matrix* vTmp = mMul(mat, v);
        mFree(result);
        result = v;
        v = vTmp;

        double resultScalar = mRow(v, 0)[0] / mRow(result, 0)[0];
        if(fabs(resultScalar - prevScalar) < 1e-7) {
            break;
        }
        prevScalar = resultScalar;
    }
    mFree(v);
    mFree(result);
    return prevScalar;
}

Matrix* mCalcEigVecSymm(Matrix* mat, double maxEig) {
    assert(mat->rows == mat->cols);
    Matrix* mat2 = mCopy(mat);
    for(int i = 0; i < mat2->rows; i++) {
        mRow(mat2, i)[i] -= maxEig;
    }

    for(int i = 0; i < mat2->rows - 1; i++) {
        double diag = mRow(mat2, i)[i];
        for(int j = i; j < mat2->cols; j++) {
            mRow(mat2, i)[j] /= diag;
        }
        for(int j = 0; j < mat2->rows; j++) {
            if(j == i)
                continue;
            double coeff = mRow(mat2, j)[i];
            for(int k = 0; k < mat->cols; k++) {
                mRow(mat2, j)[k] -= coeff * mRow(mat2, i)[k];
            }
        }
    }

    Matrix* result = mAlloc(mat2->rows, 1);
    for(int i = 0; i < mat2->rows - 1; i++) {
        mRow(result, i)[0] = -mRow(mat2, i)[mat2->cols - 1];
    }
    mRow(result, mat2->rows - 1)[0] = 1.;

    mFree(mat2);
    return result;
}

Matrix* mCalcEquivalentMinorMatrix(Matrix* mat, Matrix* x) {
    Matrix* result = mCreate(mat->rows - 1, mat->cols - 1);
    Matrix* hermitian = mIdentity(mat->rows);
    Matrix* hermitianInv = mIdentity(mat->rows);

    double x0Inv = mRow(hermitian, 0)[0] = 1. / mRow(x, 0)[0];
    for(int i = 1; i < hermitian->rows; i++) {
        mRow(hermitian, i)[0] = -mRow(x, i)[0] * x0Inv;
    }

    for(int i = 0; i < hermitianInv->rows; i++) {
        mRow(hermitianInv, i)[0] = -mRow(x, i)[0];
    }
    mRow(hermitianInv, 0)[0] *= -1.f;

    Matrix* m1 = mMul(hermitian, mat);
    Matrix* m2 = mMul(m1, hermitianInv);

    for(int i = 1; i < m2->rows; i++) {
        for(int j = 1; j < m2->cols; j++) {
            mRow(result, i - 1)[j - 1] = mRow(m2, i)[j];
        }
    }
    mFree(m2);
    mFree(m1);
    mFree(hermitian);
    mFree(hermitianInv);
    return result;
}

// doesn't work for all matrices, now only 3x3
void mSVD(Matrix* mat, Matrix** u, Matrix** s, Matrix** v) {
    *s = mCreate(mat->rows, mat->cols);
    *v = mCreate(mat->cols, mat->cols);

    Matrix* transpose = mCreateTranspose(mat);
    Matrix* ata = mMul(transpose, mat);
    mFree(transpose);

    Matrix* reduced = mCopy(ata);
    for(int i = 0; i < mat->cols; i++) {
        double maxEig = mCalcMaxEigenSymm(reduced);
        mRow(*s, i)[i] = maxEig;

        Matrix* eigVec = mCalcEigVecSymm(reduced, maxEig);

        Matrix* tmp = mCalcEquivalentMinorMatrix(reduced, eigVec);

        mFree(eigVec);
        mFree(reduced);
        reduced = tmp;
        tmp = NULL;
    }
    mFree(reduced);

    for(int i = 0; i < mat->cols; i++) {
        Matrix* eigVec = mCalcEigVecSymm(ata, mRow(*s, i)[i]);
        mNormalize(eigVec);

        for(int j = 0; j < eigVec->rows; j++) {
            mRow(*v, j)[i] = mRow(eigVec, j)[0];
        }
        mRow(*s, i)[i] = sqrt(fabs(mRow(*s, i)[i]));
        mFree(eigVec);
    }

    Matrix* matv = mMul(mat, *v);
    Matrix* sInv = mCopy(*s);
    for(int i = 0; i < sInv->cols; i++) {
        mRow(sInv, i)[i] = 1. / mRow(*s, i)[i];
    }
    *u = mMul(matv, sInv);
    mRow(*u, 0)[2] = (mRow(*u, 1)[0] * mRow(*u, 2)[1] - mRow(*u, 2)[0] * mRow(*u, 1)[1]);
    mRow(*u, 1)[2] = -(mRow(*u, 0)[0] * mRow(*u, 2)[1] - mRow(*u, 2)[0] * mRow(*u, 0)[1]);
    mRow(*u, 2)[2] = -(mRow(*u, 0)[1] * mRow(*u, 1)[0] - mRow(*u, 0)[0] * mRow(*u, 1)[1]);
    mFree(matv);
    mFree(sInv);

    mFree(ata);
}

Matrix* mCalcEigVecSymmIterative(Matrix* mat) {
    Matrix* v = mCreate(mat->rows, 1);
    for(int i = 0; i < mat->cols; i++) {
        mRow(v, i)[0] = 1.;
    }
    mNormalize(v);
    Matrix* v1 = NULL;
    for(;;) {
        v1 = mMul(mat, v);
        mNormalize(v1);
        if(mIsClose(v1, v)) {
            break;
        }
        mFree(v);
        v = v1;

        mNormalize(v);
    }
    mFree(v1);

    return v;
}

Matrix* mCalcEigensSymmetric3Cardano(Matrix* mat) {
    assert(mat->rows == mat->cols && mat->rows == 3);
    
    Matrix* eigens = mAlloc(3, 1);

    double p1 = mRow(mat, 0)[1] * mRow(mat, 0)[1]
        + mRow(mat, 0)[2] * mRow(mat, 0)[2]
        + mRow(mat, 1)[2] * mRow(mat, 1)[2];
    if( p1 < 1e-7) {
        mRow(eigens, 0)[0] = mRow(mat, 0)[0];
        mRow(eigens, 1)[1] = mRow(mat, 1)[1];
        mRow(eigens, 2)[2] = mRow(mat, 2)[2];
        return eigens;
    }
    double q = mTrace(mat) / 3.;

    double x1 = mRow(mat, 0)[0] - q;
    double x2 = mRow(mat, 1)[1] - q;
    double x3 = mRow(mat, 2)[2] - q;

    double p2 = x1*x1 + x2*x2 + x3*x3 + 2 * p1; 
    double p = sqrt(p2 / 6.);

    Matrix* B = mCopy(mat);
    for(int i = 0; i < 3; i++) {
        mRow(B, i)[i] -= q;
    }
    mScale(B, 1. / p);
    double r = mDet(B) / 2.;
    mFree(B);

    double phi = acos(r) / 3.;
    if(r <= -1.) {
        phi = M_PI / 3.;
    } else if(r >= 1.) {
        phi = 0.;
    }

    double eig1 = mRow(eigens, 0)[0] = q + 2. * p * cos(phi);
    double eig3 = mRow(eigens, 2)[0] = q + 2. * p * cos(phi + (2. * M_PI / 3.));
    mRow(eigens, 1)[0] = 3 * q - eig1 - eig3;

    if(mRow(eigens, 0)[0] < mRow(eigens, 1)[0]) {
        swap(&mRow(eigens, 0)[0], &mRow(eigens, 1)[1]);
    }
    if(mRow(eigens, 1)[0] < mRow(eigens, 2)[0]) {
        swap(&mRow(eigens, 1)[0], &mRow(eigens, 2)[0]);
    }
    if(mRow(eigens, 0)[0] < mRow(eigens, 1)[0]) {
        swap(&mRow(eigens, 0)[0], &mRow(eigens, 1)[1]);
    }

    return eigens;
}

void mSVD3(Matrix* mat, Matrix** u, Matrix** s, Matrix** v) {
    assert(mat->rows == mat->cols && mat->rows == 3);

    Matrix* at = mCreateTranspose(mat);

    Matrix* ata = mMul(at, mat);
    Matrix* aat = mMul(mat, at);
}