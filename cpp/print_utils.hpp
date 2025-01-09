#ifndef PRINT_UTILS_HPP
#define PRINT_UTILS_HPP

#include <vector>
#include <cstdio>

// Print matrix in numpy style
inline void print_matrix(const double* mat, int rows, int cols, const char* label = nullptr, const char* fmt = "%g") {
    if(label) printf("%s\n", label);
    printf("[");
    for(int i = 0; i < rows; i++) {
        printf("[");
        for(int j = 0; j < cols; j++) {
            printf(fmt, mat[i * cols + j]);
            if(j < cols-1) printf(" ");
        }
        printf("]");
        if(i < rows-1) printf("\n ");
    }
    printf("]\n");
}

// Print vector as nx1 matrix
inline void print_vector(const double* vec, int size, const char* label = nullptr) {
    print_matrix(vec, size, 1, label);
}

// Print vector of vectors as matrix
inline void print_vector(const std::vector<std::vector<int>>& vec, const char* label = nullptr) {
    if(label) printf("%s", label);
    printf("[");
    for(int i = 0; i < vec.size(); i++) {
        printf("[");
        for(int j = 0; j < vec[i].size(); j++) {
            printf("%d", vec[i][j]);
            if (j < vec[i].size()-1) printf(" ");
        }
        printf("]");
        if (i < vec.size()-1) printf("\n ");
    }
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
