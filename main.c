#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrix.h"

int main(int argc, char** argv) {
    Matrix* a = mCreate(3, 3);
    mRow(a, 0)[0] = 1.166666666666666;
    mRow(a, 0)[1] = -0.3333333333333333;
    mRow(a, 0)[2] = -0.8333333333333334;

    mRow(a, 1)[0] = -0.3333333333333333;
    mRow(a, 1)[1] = 0.6666666666666666;
    mRow(a, 1)[2] = -0.3333333333333333;

    mRow(a, 2)[0] = -0.8333333333333334;
    mRow(a, 2)[1] = -0.3333333333333333;
    mRow(a, 2)[2] = 1.1666666666666667;

    puts("A:");
    mPrint(a);
    Matrix* u = NULL;
    Matrix* s = NULL;
    Matrix* v = NULL;
    mSVD(a, &u, &s, &v);

    puts("Result U:");
    mPrint(u);
    puts("Result S:");
    mPrint(s);
    puts("Result V:");
    mPrint(v);
    

    // Matrix* m1 = mMul(u, s);
    // Matrix* m2 = mMul(m1, vt);
    // mPrint(m2);

    Matrix* uInv = mCopy(s);

    for(int i = 0; i < s->rows; i++) {
        mRow(uInv, i)[i] = 1. / mRow(uInv, i)[i];
    }
    puts("INV");
    mPrint(uInv);

    Matrix *ut = mCreateTranspose(u);
    Matrix *vt = mCreateTranspose(v);

    Matrix* b = mCreate(3, 1);
    mRow(b, 0)[0] = -1;
    mRow(b, 1)[0] = 0;
    mRow(b, 2)[0] = 1;

    Matrix* m1 = mMul(v, uInv);
    Matrix* m2 = mMul(m1, ut);
    Matrix* m3 = mMul(m2, b);

    puts("X:");
    mPrint(m3);

    Matrix* r = mMul(a, m3);
    puts("A*X");
    mPrint(r);


    return 0;
}

