#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "blas.h"


// UTILS FUNCTIONS    
double* initRandomVector(int n, int min, int max){
    double* x = calloc(n, sizeof(double));
    for (int i = 0; i < n; i++){
        x[i] = (rand() % (min - max + 1)) + min;
    }
    return x;
}

double* initRandomMatrix(double* A, int row,int col, int min, int max){
    if(A == NULL)
        A = calloc(row*col, sizeof(double));
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
            A[i*col+j] = (rand() % (min - max + 1)) + min;
        }
    }
    return A;
}

void displayVector(double* x, int n){
    printf("\n");
    for (int i = 0; i < n; i++){
            printf("%f ", x[i]);
    }
    printf("\n");
}

void displayMatrix(int row, int col, double* matrix){
    printf("_______________________________\n");
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
            printf("%0.2f  ", matrix[i*col+j]);
        }
        printf("\n");
    }
    printf("_______________________________\n");
}

double* transposeMatrix(int row, int col, double* A){
    double* AT = calloc(row*col, sizeof(double));
    for(int i=0; i<row; i++) {
        for(int j=0; j<col; j++) {
            AT[j*row+i] = A[i*col+j];
        }
    }
    return AT;
}

void Vector_Vector_Addition(double* x, double* y, int n){
    for (int i = 0; i < n; i++){
        y[i] = x[i] + y[i];
    }
}

// Produit scalaire vecteur
void Vector_Scalar_Product(double* x, double alpha, int n) {
    for (int i = 0; i < n; i++){
        x[i] *=  alpha;
    }
}


// Produit scalaire 2 vecteurs
double DotProduct(double* x, double* y, int n) {
    double result = 0;
    for (int i = 0; i < n; i++){
        result += x[i] * y[i];
    }
    return result;
}

// Produit scalaire 2 vecteurs
double DotProduct_parallel(double* x, double* y, int n) {
    double result = 0;
    // int counter[4]= {0};
    #pragma omp parallel
    {
        // printf("Thread ID == %d\n", omp_get_thread_num());
        #pragma omp for schedule(static) reduction(+:result)
            for (int i = 0; i < n; i++){
                result += x[i] * y[i];
                // printf("ThreadID : %d, i == %d\n", omp_get_thread_num(), i);
                // counter[omp_get_thread_num()] += i;
            }
    }
    // for (int i = 0; i < 4; i++){
    //     printf("Counters[%d] == %d\n", i, counter[i]);
    // }
    return result;
}

void Matrix_Vector_Product(double* A, double* v, int row, int col, double* Av){
    for(int i=0; i<row; i++) {
        for(int j=0; j<col; j++) {
            Av[i] += A[i*col+j] * v[j];
        }
    }
}

void Matrix_Vector_Product_parralel(double* A, double* v, int row, int col, double* Av){
    omp_set_nested(2);
    double Avi;
    #pragma omp parallel for schedule(static)
    for(int i=0; i<row; i++){
        Avi = Av[i];
        #pragma omp parallel for schedule(static) reduction(+:Avi)
        for(int j=0; j<col; j++) {
            Avi += A[i*col+j] * v[j];
        }
        Av[i] = Avi;
    }
}

double* Matrix_Matrix_Product(double* A, double* B, int rowA, int colA, int colB, double* AB){
    if(AB == NULL)
        AB = calloc(rowA*colB, sizeof(double));
    for(int i=0; i<rowA; i++) {
        for(int k=0; k<colB; k++) {
            for(int j=0; j<colA; j++) {
                AB[i*colB+k] += A[i*colA+j] * B[j*colB+k];
            }
        }
    }
    return AB;
}


double* Matrix_Matrix_Subsctraction(double* A, double* B, int row, int col, double* AB){
    if(AB == NULL)
        AB = calloc(row*col, sizeof(double));
    for(int i=0; i<row; i++) {
        for(int k=0; k<col; k++) {
            AB[i*col+k] = A[i*col+k] - B[i*col+k];
        }
    }
    return AB;
}

// Norme2 d'un vector
double Norme(double* x, int n){
    return sqrt(DotProduct(x, x, n));
}

// La norme Frobenius d’une matrice
double NormeFrobenius(int nb_ligne, int nb_colonne, double* matrice){
    double result = 0;
    for(int i=0; i<nb_ligne; i++) {
        for(int j=0; j<nb_colonne; j++) {
            result += pow(matrice[i*nb_colonne+j], 2);
        }
    }
    return sqrt(result);
} 
