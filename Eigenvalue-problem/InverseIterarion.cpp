#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "Norms.h"
#include "Gauss.h"

extern double EPS;

void InverseIterations(double** A_orig, int n, double* eigenvalues, double** eigenvectors) {
    std::ofstream out("InverseIterations_eigenvectors.txt");

    double* x = new double[n + 1];
    double* y = new double[n + 1];

    //M = A - lambda * E
    double** M = new double* [n + 1];
    for (int i = 1; i < n + 1; ++i) {
        M[i] = new double[n + 1];
    }

    for (int idx = 1; idx < n + 1; ++idx) {
        double lambda = eigenvalues[idx];

        for (int i = 1; i < n + 1; ++i) {
            for (int j = 1; j < n + 1; ++j) {
                M[i][j] = A_orig[i][j];
            }
            M[i][i] -= lambda;
        }

        //Initial approximation
        for (int i = 1; i < n + 1; ++i) {
            x[i] = 1.0;
        }
        double norm_x = vector_norm_2(x, n);
        for (int i = 1; i < n + 1; ++i) {
            x[i] /= norm_x;
        }

        int step = 0;
        double diff = 1.0;

        while (diff > EPS && step < 100) {
            step++;

            //Solving (A - lambda*E) * y^(k+1) = x^(k)
            Gauss_method(M, y, x, n);

            //Normalize the vector y
            double norm_y = vector_norm_2(y, n);

            //Checking convergence
            diff = 0.0;
            for (int i = 1; i <= n; ++i) {
                double next_x = y[i] / norm_y;
                //Checking the minimum distance considering the sign ambiguity of eigenvectors
                diff += std::min(std::pow(x[i] - next_x, 2), std::pow(x[i] + next_x, 2));
            }
            diff = std::sqrt(diff);

            for (int i = 1; i <= n; ++i) {
                x[i] = y[i] / norm_y;
            }
        }

        for (int i = 1; i <= n; ++i) {
            eigenvectors[i][idx] = x[i];
        }

        out << "Eigenvalue = " << eigenvalues[idx] << ", eigenvector (step = " << step << ")\n";
        for (int i = 1; i <= n; ++i) {
            out << x[i] << " ";
        }
        out << "\n\n";
    }
    out.close();

    for (int i = 1; i <= n; ++i) {
        delete[] M[i];
    }
    delete[] M;
    delete[] x;
    delete[] y;
}

void VerifyEigenvectors(double** A_orig, int n, double* eigenvalues, double** eigenvectors) {
    double* Ax = new double[n + 1];

    for (int idx = 1; idx <= n; ++idx) {
        double lambda = eigenvalues[idx];

        //Ax = A * x
        for (int i = 1; i <= n; ++i) {
            Ax[i] = 0.0;
            for (int j = 1; j <= n; ++j) {
                Ax[i] += A_orig[i][j] * eigenvectors[j][idx];
            }
        }

        //R_i = Ax_i - lambda * x_i
        double residual_norm_sq = 0.0;
        for (int i = 1; i <= n; ++i) {
            double r = Ax[i] - lambda * eigenvectors[i][idx];
            residual_norm_sq += r * r;
        }
        double residual_norm = std::sqrt(residual_norm_sq);

        std::cout << "Eigenvalue = " << lambda << "\n";
        std::cout << "Residual norm = " << residual_norm <<"\n\n";
    }
    delete[] Ax;
}

void RayleighIteration(double** A, int n, double* initial_eigenvalues, double* eigenvalues, double** eigenvectors) {
    std::ofstream out("Rayleigh_eigensystem.txt");

    double* x = new double[n + 1];
    double* y = new double[n + 1];
    double* Ax = new double[n + 1];

    double** M = new double* [n + 1];
    for (int i = 1; i < n + 1; ++i) {
        M[i] = new double[n + 1];
    }

    int max_iterations = 100;

    for (int idx = 1; idx <= n; ++idx) {
        double lambda = initial_eigenvalues[idx];
        double lambda_old = lambda;

        for (int i = 1; i <= n; ++i) {
            x[i] = 1.0;
        }
        double x0_norm = vector_norm_2(x, n);
        for (int i = 1; i <= n; ++i) {
            x[i] /= x0_norm;
        }

        int step = 0;

        for (int k = 1; k <= max_iterations; ++k) {
            //M = A - lambda * E
            for (int i = 1; i <= n; ++i) {
                for (int j = 1; j <= n; ++j) {
                    M[i][j] = A[i][j];
                }
                M[i][i] -= lambda;
            }

            Gauss_method(M, y, x, n);

            double norm_y = vector_norm_2(y, n);
            if (norm_y < 1e-15) {
                break;
            }
            for (int i = 1; i <= n; ++i) {
                x[i] = y[i] / norm_y;
            }

            for (int i = 1; i <= n; ++i) {
                Ax[i] = 0.0;
                for (int j = 1; j <= n; ++j) {
                    Ax[i] += A[i][j] * x[j];
                }
            }

            //lambda = (Ax, x) / (x, x)
            double sum = 0.0;
            for (int i = 1; i <= n; ++i) {
                sum += Ax[i] * x[i];
            }
            lambda = sum;

            //Lambda convergence
            if (k > 1 && std::fabs(lambda - lambda_old) < EPS) {
                step = k;
                break;
            }
            lambda_old = lambda;
            step = k;
        }

        eigenvalues[idx] = lambda;
        for (int i = 1; i <= n; ++i) {
            eigenvectors[i][idx] = x[i];
        }

        out << "Eigenvalue = " << lambda << " (steps = " << step << ")\n";
        out << "Eigenvector:\n";
        for (int i = 1; i <= n; ++i) {
            out << x[i] << " ";
        }
        out << "\n\n";
    }

    out.close();

    for (int i = 1; i <= n; ++i) {
        delete[] M[i];
    }
    delete[] M;
    delete[] x;
    delete[] y;
    delete[] Ax;
}