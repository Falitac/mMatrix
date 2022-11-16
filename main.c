#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"


float nextNums(Matrix* m, int i, int j) {
    return (i==j) * 4 + i + 1 + j * j;
}

int main(int argc, char** argv) {
    Matrix* m = mCreate(3, 3);
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
            mRow(m, i)[j] = i * 3 + j;
    mPrint(m);
    Matrix* n = mCreateTranspose(m);
    mPrint(n);
    Matrix* out = mMul(m, n);
    mScale(out, 0.25);
    mAdd(out, out);
    mGenerate(out, nextNums);

    mPrint(out);
    printf("Det: %.3f\n", mDet(out));
    mRedRows(out);
    mPrint(out);
    printf("Fact. diagonal: %.3f\n", mFactorDiagonal(out));
    mTranspose(out);
    mPrint(out);
    mRedRows(out);
    mPrint(out);

    mFree(out);
    mFree(n);
    mFree(m);
    return 0;
}
