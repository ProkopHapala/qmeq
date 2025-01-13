#ifndef GAUSS_SOLVER_HPP
#define GAUSS_SOLVER_HPP

#include <cmath>

// class GaussSolver {
// public:
//     // Solve the system Ax = b using Gaussian elimination with partial pivoting
//     // Returns the solution in the x array
//     static void solve(double* A, double* b, double* x, int n) {
//         double* aug_matrix = new double[n * (n + 1)];  // Augmented matrix [A|b]
        
//         // Create augmented matrix
//         for(int i = 0; i < n; i++) {
//             for(int j = 0; j < n; j++) {
//                 aug_matrix[i * (n + 1) + j] = A[i * n + j];
//             }
//             aug_matrix[i * (n + 1) + n] = b[i];
//         }

//         // Gaussian elimination with partial pivoting
//         for(int k = 0; k < n - 1; k++) {
//             // Find pivot
//             int pivot_row = k;
//             double pivot_val = fabs(aug_matrix[k * (n + 1) + k]);
//             for(int i = k + 1; i < n; i++) {
//                 double val = fabs(aug_matrix[i * (n + 1) + k]);
//                 if(val > pivot_val) {
//                     pivot_val = val;
//                     pivot_row = i;
//                 }
//             }

//             // Swap rows if necessary
//             if(pivot_row != k) {
//                 for(int j = k; j <= n; j++) {
//                     double temp = aug_matrix[k * (n + 1) + j];
//                     aug_matrix[k * (n + 1) + j] = aug_matrix[pivot_row * (n + 1) + j];
//                     aug_matrix[pivot_row * (n + 1) + j] = temp;
//                 }
//             }

//             // Eliminate column
//             for(int i = k + 1; i < n; i++) {
//                 double factor = aug_matrix[i * (n + 1) + k] / aug_matrix[k * (n + 1) + k];
//                 for(int j = k; j <= n; j++) {
//                     aug_matrix[i * (n + 1) + j] -= factor * aug_matrix[k * (n + 1) + j];
//                 }
//             }
//         }

//         // Back substitution
//         for(int i = n - 1; i >= 0; i--) {
//             x[i] = aug_matrix[i * (n + 1) + n];
//             for(int j = i + 1; j < n; j++) {
//                 x[i] -= aug_matrix[i * (n + 1) + j] * x[j];
//             }
//             x[i] /= aug_matrix[i * (n + 1) + i];
//         }

//         delete[] aug_matrix;
//     }
// };


void GaussElimination( int n, double ** A, double * c, int * index ) {

	// Initialize the index
	for (int i=0; i<n; ++i) index[i] = i;

	// Find the rescaling factors, one from each row
	for (int i=0; i<n; ++i) {
	  double c1 = 0;
	  for (int j=0; j<n; ++j) {
		double c0 = fabs( A[i][j] );
		if (c0 > c1) c1 = c0;
	  }
	  c[i] = c1;
	}

	// Search the pivoting element from each column
	int k = 0;
	for (int j=0; j<n-1; ++j) {
	  double pi1 = 0;
	  for (int i=j; i<n; ++i) {
		double pi0 = fabs( A[ index[i] ][j] );
		pi0 /= c[ index[i] ];
		if (pi0 > pi1) {
		  pi1 = pi0;
		  k = i;
		}
	  }

	  // Interchange rows according to the pivoting order
	  int itmp = index[j];
	  index[j] = index[k];
	  index[k] = itmp;
	  for (int i=j+1; i<n; ++i) {
		double pj = A[ index[i] ][ j ]/A[ index[j] ][j];

	   // Record pivoting ratios below the diagonal
		A[ index[i] ][j] = pj;

	   // Modify other elements accordingly
		for (int l=j+1; l<n; ++l)
		  A[ index[i] ][l] -= pj*A[ index[j] ][l];
	  }
	}
}

void linSolve_gauss( int n, double ** A, double * b, int * index, double * x ) {

	// Transform the matrix into an upper triangle
	GaussElimination( n, A, x, index);

	// Update the array b[i] with the ratios stored
	for(int i=0; i<n-1; ++i) {
		for(int j =i+1; j<n; ++j) {
			b[index[j]] -= A[index[j]][i]*b[index[i]];
		}
	}

	// Perform backward substitutions
	x[n-1] = b[index[n-1]]/A[index[n-1]][n-1];
	for (int i=n-2; i>=0; --i) {
		x[i] = b[index[i]];
		for (int j=i+1; j<n; ++j) {
			x[i] -= A[index[i]][j]*x[j];
		}
		x[i] /= A[index[i]][i];
	}
}

void linSolve_gauss( int n, double* A, double* b, double * x ) {
    int index[n];
    double* A_[n];
    for(int i = 0; i < n; i++) { A_[i] = &A[i * n]; }
    linSolve_gauss( n, A_, b, index, x );
}

#endif // GAUSS_SOLVER_HPP
