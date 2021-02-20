#include <iostream>
#include "FORSYTHE.h"
#include "MATRIX.H"

// Оценка числа обусловленности
N_TYPE cond;

// Создание пустой матрицы (ndim,ndim)
N_TYPE **makeEmptyMatrix() {
    auto **Matrix = new N_TYPE *[ndim];
    for (int i = 0; i < ndim ; ++i) Matrix[i] = new N_TYPE [ndim];
    return Matrix;
}

// Умножение двух матриц
N_TYPE **multiply(N_TYPE **A, N_TYPE **invA) {
    auto **RES = makeEmptyMatrix();
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) RES[i][j] = 0.0;
    }
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
            for (int k = 0; k < ndim; ++k) RES[i][j] += A[i][k] * invA[k][j];
        }
    }
    return RES;
}

// Вычитание одной матрицы из другой
N_TYPE **minus(N_TYPE **A, N_TYPE **B) {
    auto **RES = makeEmptyMatrix();
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) RES[i][j] = A[i][j] - B[i][j];
    }
    return RES;
}

// Создание матрицы A
N_TYPE **makeAMatrix(double x) {
    auto **A = makeEmptyMatrix();
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
            if (i == j) A[i][j] = 1 + x;
            else A[i][j] = j + 1.0;
        }
    }
    A[0][0] = 1.0;
    return A;
}

// Создание единичной матрицы
N_TYPE **makeEMatrix() {
    auto **E = makeEmptyMatrix();
    for (int i = 0; i < ndim ; ++i) {
        for (int j = 0; j < ndim; ++j) {
            if (i == j) E[i][j] = 1.0;
            else E[i][j] = 0.0;
        }
    }
    return E;
}

// Создание обратной матрицы A
N_TYPE **makeInverseA(N_TYPE **A) {
    N_TYPE Work[ndim], subMatrix[ndim][ndim];
    int IPVT[ndim];
    auto **invA = makeEmptyMatrix();
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) subMatrix[i][j] = A[i][j];
    }

    decomp(ndim, subMatrix, &cond, IPVT, Work);

    for (int i = 0; i < ndim; ++i) {
        N_TYPE B[ndim]{};
        B[i] = 1.0;
        solve(ndim, subMatrix, B, IPVT);
        for (int j = 0; j < ndim; ++j) invA[j][i] = B[j];
    }

    return invA;
}

// Создание матрицы R
N_TYPE **makeRMatrix(N_TYPE **A, N_TYPE **invA, N_TYPE **E) {
    N_TYPE **R;
    N_TYPE **AinvA = multiply(A, invA);
    R = minus(AinvA, E);
    return R;
}

// Вычисление нормы матрицы
double norm(N_TYPE **A) {
    double norm = 0;
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) norm += sqrt(pow(A[i][j], 2));
    }
    return norm;
}


// Вывод матрицы в консоль
void print(N_TYPE **matrix)  {
    for (int i = 0; i < ndim ; ++i) {
        for (int j = 0; j < ndim; ++j) printf("%.15lf      ", matrix[i][j]);
        printf("\n");
    }
}

int main() {
    double x = 1.1;
    N_TYPE **A = makeAMatrix(x);
    N_TYPE **invA = makeInverseA(A);
    N_TYPE **E = makeEMatrix();
    N_TYPE **R = makeRMatrix(A, invA, E);
    printf("A: \n");
    print(A);
    printf("\nA^-1: \n");
    print(invA);
    printf("\nR: \n");
    print(R);
    double norm1 = norm(R);
    printf("\n   x = %lf\n", x);
    printf("cond = %.15lf\n", cond);
    printf("norm = %.15lf", norm1);
}
