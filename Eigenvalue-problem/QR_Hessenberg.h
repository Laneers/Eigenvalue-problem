#pragma once

void Hessenberg(double** A, int n, bool isSymmetric = false);
void QRIterationShift(double** A, int m, int n);
void EigenvaluesQR(double** A, int n, double* eigenvalues);
void verifyEigenvalues(double** A_original, int n, double* eigenvalues);