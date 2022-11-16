#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

int main(int argc, char** argv) {
    Matrix* m = mCreate(3, 3);
    for(int i = 0; i < 3; i++) 
        for(int j = 0; j < 3; j++) 
            mRow(m, i)[j] = i * 3 + j;
    mPrint(m);
    Matrix* n = mCreateTranspose(m);
    mPrint(n);
    Matrix* out = mMultiply(m, n);
    mScale(out, 0.25);
    mPrint(out);

    mFree(out);
    mFree(n);
    mFree(m);
    return 0;
}
