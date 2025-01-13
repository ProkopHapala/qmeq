#ifndef PRINT_UTILS_HPP
#define PRINT_UTILS_HPP

#include <vector>
#include <cstdio>

inline void print_vector(const double* vec, int n, const char* fmt = "%g ") {
    printf("[ "); 
    for(int j = 0; j < n; j++) { printf(fmt, vec[j]); }
    printf("]\n"); 
}

inline void print_vector(const int* vec, int n, const char* fmt = "%i ") {
    printf("[ "); 
    for(int j = 0; j < n; j++) { printf(fmt, vec[j]); }
    printf("],"); 
}
inline void print_vector( const std::vector<int>& vec, const char* fmt = "%i ") { print_vector( vec.data(), vec.size(), fmt); }

// Print matrix in numpy style
inline void print_matrix(const double* mat, int rows, int cols, const char* fmt = "%g ") {
    printf("[");
    for(int i = 0; i < rows; i++) { print_vector( mat+i*cols, cols, fmt); }
    printf("]\n");
}

// Print vector as nx1 matrix
// inline void print_vector(const double* vec, int size, const char* label = nullptr) {
//     print_matrix(vec, size, 1, label);
// }

// Print vector of vectors as matrix
inline void print_vector_of_vectors(const std::vector<std::vector<int>>& vec, const char* fmt = "%i " ) {
    //if(label) printf("%s", label);
    printf("[");
    for(int i = 0; i < vec.size(); i++) { print_vector( vec[i] ); }
    printf("]\n");
}

// Print 3D array as a sequence of matrices
inline void print_3d_array(const double* arr, int dim1, int dim2, int dim3, const char* label = "Matrix") {
    //if(label) printf("%s", label);
    for(int i = 0; i < dim1; i++) {
        printf("\n%s %i:\n", label, i);
        print_matrix(&arr[i * dim2 * dim3], dim2, dim3);
    }
}

#endif // PRINT_UTILS_HPP
