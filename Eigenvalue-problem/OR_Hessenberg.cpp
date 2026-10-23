#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

extern double EPS;

void Hessenberg(double** A, int n, bool isSymmetric = false) {
    std::ofstream out("Hess_sequence.txt");
    for (int j = 1; j < n - 1; ++j) {
        int k = j + 1;

        for (int l = k + 1; l < n + 1; ++l) {
            if (std::fabs(A[l][j]) < 1e-15) continue;

            double denominator = std::sqrt(A[k][j] * A[k][j] + A[l][j] * A[l][j]);
            double alpha = A[k][j] / denominator;
            double beta = A[l][j] / denominator;

            //Evaluating only lower triangle
            if (isSymmetric) {
                //Left roatations
                for (int col = 1; col < k + 1; ++col) {
                    double tmp_k = alpha * A[k][col] + beta * A[l][col];
                    double tmp_l = -beta * A[k][col] + alpha * A[l][col];
                    A[k][col] = tmp_k;
                    A[l][col] = tmp_l;
                }

                double akk = A[k][k], akl = A[l][k], all = A[l][l];

                A[k][k] = alpha * alpha * akk + 2.0 * alpha * beta * akl + beta * beta * all;
                A[l][l] = beta * beta * akk - 2.0 * alpha * beta * akl + alpha * alpha * all;
                A[l][k] = -beta * alpha * akk + (alpha * alpha - beta * beta) * akl + alpha * beta * all;

                //Right rotation
                for (int row = l + 1; row < n + 1; ++row) {
                    double tmp_k = alpha * A[row][k] + beta * A[row][l];
                    double tmp_l = -beta * A[row][k] + alpha * A[row][l];
                    A[row][k] = tmp_k;
                    A[row][l] = tmp_l;
                }
                A[l][j] = 0.0;

                //Restoring upper half
                for (int row = 1; row < n + 1; ++row) {
                    for (int col = row + 1; col < n + 1; ++col) {
                        A[row][col] = A[col][row];
                    }
                }
            }
            else {
                //Left rotation
                for (int col = 1; col < n + 1; ++col) {
                    double tmp_k = alpha * A[k][col] + beta * A[l][col];
                    double tmp_l = -beta * A[k][col] + alpha * A[l][col];
                    A[k][col] = tmp_k;
                    A[l][col] = tmp_l;
                }

                //Right rotation
                for (int row = 1; row < n + 1; ++row) {
                    double tmp_k = alpha * A[row][k] + beta * A[row][l];
                    double tmp_l = -beta * A[row][k] + alpha * A[row][l];
                    A[row][k] = tmp_k;
                    A[row][l] = tmp_l;
                }

                A[l][j] = 0.0;
            }
            std::string stepInfo = "Resetting an element A[" + std::to_string(l) + "][" + std::to_string(j) + "]:\n";
            out << stepInfo;
            for (int i = 1; i < n + 1; i++) {
                for (int j = 1; j < n + 1; j++) {
                    out << A[i][j] << " ";
                }
                out << "\n";
            }
        }
    }
    out.close();
}

void QRIterationShift(double** A, int m, int n) {
    double sigma = A[m - 1][m - 1];

    for (int i = 1; i < m + 1; ++i) {
        A[i][i] -= sigma;
    }

    double* c = new double[m - 1];
    double* s = new double[m - 1];

    for (int i = 1; i < m; ++i) {
        int k = i;
        int l = i + 1;
        if (std::fabs(A[l][i]) > 1e-15) {
            double denominator = std::sqrt(A[k][i] * A[k][i] + A[l][i] * A[l][i]);
            c[i] = A[k][i] / denominator;
            s[i] = A[l][i] / denominator;

            //Left
            for (int col = i; col < m + 1; ++col) {
                double tmp_k = c[i] * A[k][col] + s[i] * A[l][col];
                double tmp_l = -s[i] * A[k][col] + c[i] * A[l][col];
                A[k][col] = tmp_k;
                A[l][col] = tmp_l;
            }
        }
        else {
            c[i] = 1.0;
            s[i] = 0.0;
        }
    }

    for (int i = 1; i < m; ++i) {
        int k = i;
        int l = i + 1;
        //Right
        for (int row = 1; row < m + 1; ++row) {
            double tmp_k = c[i] * A[row][k] + s[i] * A[row][l];
            double tmp_l = -s[i] * A[row][k] + c[i] * A[row][l];
            A[row][k] = tmp_k;
            A[row][l] = tmp_l;
        }
    }

    //Back shift
    for (int i = 1; i < m + 1; ++i) {
        A[i][i] += sigma;
    }

    delete[] c;
    delete[] s;
}

void EigenvaluesQR(double** A, int n, double* eigenvalues) {
    std::ofstream out("QR_eigenvalues.txt");

    double** H = new double* [n];
    for (int i = 0; i < n; ++i) {
        H[i] = new double[n];
        for (int j = 0; j < n; ++j) H[i][j] = A[i][j];
    }

    int m = n;
    int step = 0;

    while (m > 1) {
        while (std::fabs(H[m - 1][m - 2]) > EPS) {
            QRIterationShift(H, m, n);
            step++;
        }
        eigenvalues[m - 1] = H[m - 1][m - 1];

        out << eigenvalues[m - 1] << " (step = " << step << ")\n";
        m--;
        step = 0;
    }

    //The last eigenvalue
    eigenvalues[0] = H[0][0];
    out << eigenvalues[0] << "\n";
    out.close();

    for (int i = 0; i < n; ++i) {
        delete[] H[i];
    }
    delete[] H;
}