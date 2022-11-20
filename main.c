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
    mSetPrintNumberFormat("%7.4f");
    Matrix* a = mLoadFromFile("../in3.txt");

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

    puts("I am here, :3");
    return 0;
}
