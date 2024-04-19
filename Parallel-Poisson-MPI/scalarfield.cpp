#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#define N 120

void printScalarField(double** phi, int numRows, int numCols) {
    std::cout << numRows << " " << numCols << std::endl;
    /* Print the values of the scalar field */
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << std::fixed << std::setprecision(6) << phi[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    /* Seed the random number generator */
    std::srand(std::time(nullptr));
    double** phi = new double*[N];
    for (int i = 0; i < N; ++i) {
        phi[i] = new double[N];
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            phi[i][j] = static_cast<double>(std::rand()) / RAND_MAX;
        }
    }

    printScalarField(phi, N, N);

    /* Deallocate memory for the scalar field */
    for (int i = 0; i < N; ++i) {
        delete[] phi[i];
    }
    delete[] phi;

    return 0;
}
