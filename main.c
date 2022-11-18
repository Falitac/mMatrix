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
    Matrix* a = mLoadFromFile("../trans.txt");
    Matrix* v = mLoadFromFile("../in2v.txt");
    mNormalize(v);


    puts("A:");
    mPrint(a);

    Matrix* at = mCreateTranspose(a);
    puts("At:");
    mPrint(at);

    Matrix* aat = mMul(a, at);

    puts("A*At:");
    mPrint(aat);
    mPrint(v);

    for(int i = 0; i < 20; i++) {
        Matrix* v1 = mMul(aat, v);
        mFree(v);
        v = v1;
        v1 = NULL;

        mNormalize(v);
    }
    puts("V:");
    mPrint(v);

    return 0;
}
