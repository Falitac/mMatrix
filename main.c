#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrix.h"


float nextNums(Matrix* m, int i, int j) {
    return (i==j) * 4 + i + 1 + j * j;
}

float randUniform(float a, float b) {
    return (float) rand() / ((float)RAND_MAX / (b - a)) + a;
}

Matrix* findEigen(Matrix* mat) {
    return NULL;
}

int main(int argc, char** argv) {
    srand(time(NULL));
    Matrix* a = mLoadFromFile("../in2.txt");
    Matrix* v = mCreate(2, 1);
    mFill(v, 1.);

    mRow(v, 0)[0] = randUniform(-1.f, 1.f);
    mRow(v, 1)[0] = randUniform(-1.f, 1.f);
    //mRow(v, 2)[0] = randUniform(-1.f, 1.f);
    mNormalize(v);

    Matrix* at = mCreateTranspose(a);

    puts("A:");
    mPrint(a);

    puts("At:");
    mPrint(at);

    mMulFirst(a, at);
    mFree(at);

    puts("A*At:");
    mPrint(a);
    mPrint(v);

    for(int i = 0; i < 20; i++) {
        Matrix* v1 = mMul(a, v);
        mFree(v);
        v = v1;
        v1 = NULL;

        mNormalize(v);
    }
    puts("V:");
    mPrint(v);

    return 0;
}
