#pragma once

void InverseIterations(double** A_orig, int n, double* eigenvalues, double** eigenvectors);
void VerifyEigenvectors(double** A_orig, int n, double* eigenvalues, double** eigenvectors);
void RayleighIteration(double** A, int n, double* initial_eigenvalues, double* eigenvalues, double** eigenvectors);