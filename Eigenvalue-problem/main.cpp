#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

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

    double** A = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        A[i] = new double[n + 1];
    }

    double val;
    for (int i = 1; i < n + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            file >> val;
            A[i][j] = val;
        }
    }
    file.close();

    bool flag_sym = true;
    for (int i = 1; i < n + 1; i++) {
        for (int j = 1; j < n + 1; j++) {
            if (A[i][j] != A[j][i]) {
                flag_sym = false;
                break;
            }
        }
    }

    //Destructors
    for (int i = 1; i <= n; i++) {
        delete[] A[i];
    }
    delete[] A;
    return 0;
}