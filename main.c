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
    Matrix* a = mLoadFromFile("../in3.txt");


    puts("A:");
    mPrint(a);

    Matrix* at = mCreateTranspose(a);
    puts("At:");
    mPrint(at);

    Matrix* ata = mMul(at, a);
    Matrix* aat = mMul(a, at);

    puts("At*A:");
    mPrint(ata);
    puts("A*At:");
    mPrint(aat);

    Matrix* eig = mCalcEigensSymmetric3(aat);
    puts("Eigens:");
    mPrint(eig);

    Matrix* first = mCopy(aat);
    Matrix* second = mIdentity(aat->rows);
    mScale(second, -mRow(eig, 0)[0]);
    mAdd(first, second);

    puts("second:");
    mPrint(second);
    puts("diff:");
    mPrint(first);

    Matrix* ve = mCalcEigensSymmetric3(first);
    puts("V1:");
    mPrint(ve);

    return 0;
}
