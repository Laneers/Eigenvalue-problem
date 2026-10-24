#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "QR_Hessenberg.h"
#include "InverseIteration.h"

double EPS = 1e-6;

int main() {
    //Reading from file
    std::string filename;
    std::cin >> filename;

    int n;
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return 0;
    }
    file >> n;

    double** A_original = new double* [n + 1];
    double** A_Hessenberg = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        A_original[i] = new double[n + 1];
        A_Hessenberg[i] = new double[n + 1];
    }

    double val;
    for (int i = 1; i < n + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            file >> val;
            A_original[i][j] = val;
            A_Hessenberg[i][j] = val;
        }
    }
    file.close();
	std::cout << "\n";

    bool isSymmetric = true;
    for (int i = 1; i < n + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            if (A_original[i][j] != A_original[j][i]) {
                isSymmetric = false;
                break;
            }
        }
    }

    std::cout << "---QR iterations---\n\n";
	double* eigenvalues = new double[n + 1];
    Hessenberg(A_Hessenberg, n, isSymmetric);
    EigenvaluesQR(A_Hessenberg, n, eigenvalues);
    verifyEigenvalues(A_original, n, eigenvalues);

	double** eigenvectors = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        eigenvectors[i] = new double[n + 1];
	}

	InverseIterations(A_original, n, eigenvalues, eigenvectors);
    VerifyEigenvectors(A_original, n, eigenvalues, eigenvectors);

    double* refined_eigenvalues = new double[n + 1];
    double** refined_eigenvectors = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        refined_eigenvectors[i] = new double[n + 1];
    }

    std::cout << "---Rayleigh iterations---\n\n";
    RayleighIteration(A_original, n, eigenvalues, refined_eigenvalues, refined_eigenvectors);
    verifyEigenvalues(A_original, n, refined_eigenvalues);
    VerifyEigenvectors(A_original, n, refined_eigenvalues, refined_eigenvectors);

    //Destructors
    for (int i = 1; i <= n; i++) {
        delete[] A_original[i];
        delete[] A_Hessenberg[i];
        delete[] eigenvectors[i];
        delete[] refined_eigenvectors[i];
    }
    delete[] A_original;
    delete[] A_Hessenberg;
    delete[] eigenvalues;
    delete[] eigenvectors;
    delete[] refined_eigenvalues;
    delete[] refined_eigenvectors;

    return 0;
}