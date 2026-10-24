#include <iostream>
#include <fstream>
#include <cmath>

#include "Norms.h"

void Gauss_method(double** M, double* x, double* y, int n) {
    double** A = new double* [n + 1];
    for (int i = 1; i <= n; ++i) {
		A[i] = new double[n + 2]; //augmented matrix
        for (int j = 1; j <= n; ++j) {
            A[i][j] = M[i][j];
        }
        A[i][n + 1] = y[i];
    }

    //Selecting the main element
    int index;
    for (int i = 1; i < n; i++) {
        double maxEl = std::fabs(A[i][i]);
        index = i;
        for (int j = i + 1; j < n + 1; j++) {
            if (std::fabs(A[j][i]) > maxEl) {
                maxEl = std::fabs(A[j][i]);
                index = j;
            }
        }
        for (int k = 1; k < n + 2; k++) {     //change rows
            double temp = A[i][k];
            A[i][k] = A[index][k];
            A[index][k] = temp;
        }

        //Forward elimination
        double tmp;
        for (int l = i + 1; l < n + 1; l++) {
            tmp = A[l][i] / A[i][i];
            for (int t = i; t < n + 2; t++) {
                A[l][t] -= tmp * A[i][t];
            }
        }
    }

    //Back substitution
    x[n] = A[n][n + 1] / A[n][n];
    for (int i = n - 1; i >= 1; i--) {
        x[i] = A[i][n + 1];
        for (int j = i + 1; j <= n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] = x[i] / A[i][i];
    }

    for (int i = 1; i <= n; i++) {
        delete[] A[i];
    }
    delete[] A;
}