#ifndef GAUSS_SOLVER_HPP
#define GAUSS_SOLVER_HPP

#include <cmath>

class GaussSolver {
public:
    // Solve the system Ax = b using Gaussian elimination with partial pivoting
    // Returns the solution in the x array
    static void solve(double* A, double* b, double* x, int n) {
        double* aug_matrix = new double[n * (n + 1)];  // Augmented matrix [A|b]
        
        // Create augmented matrix
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                aug_matrix[i * (n + 1) + j] = A[i * n + j];
            }
            aug_matrix[i * (n + 1) + n] = b[i];
        }

        // Gaussian elimination with partial pivoting
        for(int k = 0; k < n - 1; k++) {
            // Find pivot
            int pivot_row = k;
            double pivot_val = fabs(aug_matrix[k * (n + 1) + k]);
            for(int i = k + 1; i < n; i++) {
                double val = fabs(aug_matrix[i * (n + 1) + k]);
                if(val > pivot_val) {
                    pivot_val = val;
                    pivot_row = i;
                }
            }

            // Swap rows if necessary
            if(pivot_row != k) {
                for(int j = k; j <= n; j++) {
                    double temp = aug_matrix[k * (n + 1) + j];
                    aug_matrix[k * (n + 1) + j] = aug_matrix[pivot_row * (n + 1) + j];
                    aug_matrix[pivot_row * (n + 1) + j] = temp;
                }
            }

            // Eliminate column
            for(int i = k + 1; i < n; i++) {
                double factor = aug_matrix[i * (n + 1) + k] / aug_matrix[k * (n + 1) + k];
                for(int j = k; j <= n; j++) {
                    aug_matrix[i * (n + 1) + j] -= factor * aug_matrix[k * (n + 1) + j];
                }
            }
        }

        // Back substitution
        for(int i = n - 1; i >= 0; i--) {
            x[i] = aug_matrix[i * (n + 1) + n];
            for(int j = i + 1; j < n; j++) {
                x[i] -= aug_matrix[i * (n + 1) + j] * x[j];
            }
            x[i] /= aug_matrix[i * (n + 1) + i];
        }

        delete[] aug_matrix;
    }
};

#endif // GAUSS_SOLVER_HPP
