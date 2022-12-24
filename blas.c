#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>
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

// Produit scalaire vecteur
void Vector_Scalar_Product_parallel(double* x, double alpha, int n) {
    #pragma omp parallel for schedule(static)
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
    #pragma omp parallel for schedule(static) reduction(+:result)
    for (int i = 0; i < n; i++){
        result += x[i] * y[i];
    }
    return result;
}

void Matrix_Vector_Product(double* A, double* v, int row, int col, double* Av, int parallel){
    if(parallel)
        Matrix_Vector_Product_parralel(A, v, row, col, Av);
    else
        Matrix_Vector_Product_sequential(A, v, row, col, Av);
}

void Matrix_Vector_Product_sequential(double* A, double* v, int row, int col, double* Av){
    for(int i=0; i<row; i++) {
        Av[i] = 0;
        for(int j=0; j<col; j++) {
            Av[i] += A[i*col+j] * v[j];
        }
    }
}

void Matrix_Vector_Product_parralel(double* A, double* v, int row, int col, double* Av){
    omp_set_nested(2);
    // double Avi;
    #pragma omp parallel for schedule(static)
    for(int i=0; i<row; i++){
        // printf("Av %d[%d]\n",omp_get_thread_num(),i);
        // Avi = Av[i];
        // #pragma omp parallel for schedule(static) reduction(+:Avi)
        double sum =  0;
        for(int j=0; j<col; j++) {
            // printf("A[%d] v %d[%d]\n",i, omp_get_thread_num(),j);
            sum += A[i*col+j] * v[j];
        }
        Av[i] = sum;
        // Av[i] = Avi;
    }
}


// Based on parralel

void blas21(double* A, double* x, double* y, double alpha, double beta, int row, int col, int nbNonZeroA, int parallel, int sparce_rep){
    if(parallel){
        if(!sparce_rep)
            blas21_parallel(A, x, y, alpha, beta, row, col);
        else
            blas21_parallel_sparce(A, x, y, alpha, beta, col, nbNonZeroA);
    }
    else{
        if(!sparce_rep)
            blas21_sequential(A, x, y, alpha, beta, row, col);
        else
            blas21_sequential_sparce(A, x, y, alpha, beta, col, nbNonZeroA);
    }
}


// x= α.Ax+βy
void blas21_sequential(double* A, double* x, double* y, double alpha, double beta, int row, int col){
    double* v = malloc(row * sizeof(double));

    for(int i=0; i<row; i++) {
        double ax = 0;
        for(int j=0; j<col; j++) {
            ax += A[i*col+j] * x[j];
        }
        v[i] = alpha * ax + beta * y[i];
    }

    memcpy(x, v, row * sizeof(double));    
    free(v);
}

// x= α.Ax+βy
void blas21_sequential_sparce(double* sparceA, double* x, double* y, double alpha, double beta, int n, int sizeSparceA){
    double* Av = calloc(n ,sizeof(double));
    // Matrix Vector Product 
    Sparce_Matrix_Vector_Product(sparceA, x, sizeSparceA, n, Av, 0);
    for(int j=0; j<n; j++) {
        x[j] = alpha * Av[j] + beta * y[j];
    }
    free(Av);
}


void blas21_parallel(double* A, double* x, double* y, double alpha, double beta, int row, int col){
    double* v = malloc(row * sizeof(double));

    #pragma omp parallel for schedule(static)
    for(int i=0; i<row; i++) {
        double ax = 0;
        for(int j=0; j<col; j++) {
            ax += A[i*col+j] * x[j];
        }
        v[i] = alpha * ax + beta * y[i];
    }

    memcpy(x, v, row * sizeof(double));    
    free(v);
}

// x= α.Ax+βy
void blas21_parallel_sparce(double* sparceA, double* x, double* y, double alpha, double beta, int n, int sizeSparceA){
    double* Av = calloc(n ,sizeof(double));
    // Matrix Vector Product 
    Sparce_Matrix_Vector_Product(sparceA, x, sizeSparceA, n, Av, 1);
    #pragma omp parallel for schedule(static)
    for(int j=0; j<n; j++) {
        x[j] = alpha * Av[j] + beta * y[j];
    }
    free(Av);
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

// Norme1 d'un vector
double Norme_One(double* x, int n){
    double result = 0;
    for (int i = 0; i < n; i++){
        result += x[i];
    }
    return result;
}

// Norme1 d'un vector
double Norme_One_parallel(double* x, int n){
    double result = 0;
    #pragma omp parallel for schedule(static) reduction(+:result)
    for (int i = 0; i < n; i++){
        result += x[i];
    }
    return result;
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




void Sparce_Matrix_Vector_Product(double* A, double* v, int sizeA, int sizeV, double* Av, int parallel){
    #pragma omp parallel if(parallel)
    {
        double* Av_private;
        if(parallel) Av_private = calloc(sizeV, sizeof(double));
        #pragma omp for schedule(static)
        for(int i=0; i<sizeA; i++) {
            int r_idx = A[i*3];
            int c_idx = A[i*3+1];
            if(parallel)
                Av_private[r_idx] += A[i*3+2] * v[c_idx];
            else
                Av[r_idx] += A[i*3+2] * v[c_idx];
        }
        if(parallel){
            #pragma omp critical
            {
                for(int i=0; i<sizeV; i++)
                    Av[i] += Av_private[i];
            }
            free(Av_private);
        }
    }
}

int count_zero_matrix(double *A, int row, int col){
	int nbNonZero = 0;
    for(int i=0; i<row; i++) {
        for(int j=0; j<col; j++) {
            if (A[i*col+j] != 0)
				nbNonZero++;
        }
    }
    return nbNonZero;
}

// Convert normal matrix representation to sparce matrix representation
double* matrix_to_sparce(double *A, int row, int col, int* nbNonZero){
	// If nbNonZero not known
	if(*nbNonZero < 0){
		*nbNonZero = count_zero_matrix(A, row, col);
	}

	double *sparceA = calloc(*nbNonZero * 3,sizeof(double));

	// Making of new matrix
	int k = 0;
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			if (A[i*col+j] != 0)
			{
				sparceA[3*k] = i;
				sparceA[3*k+1] = j;
				sparceA[3*k+2] = A[i*col+j];
				k++;
			}

	return sparceA;
}


// Convert sparce matrix representation to normal matrix representation
double* sparce_to_matrix(double *sparceA, double *A, int row, int col, int nbNonZero){
	int r_idx, c_idx;

	// Init matrix 
	if(A == NULL)
        A = calloc(row*col, sizeof(double));

	// Making of new matrix
	for (int i=0; i<nbNonZero; i++)
	{
		r_idx = sparceA[i*3];
		c_idx = sparceA[i*3+1];
		A[r_idx*col+c_idx] = sparceA[3*i+2];
	}

	return A;
}

void display_sparce_matrix(double *sparceA, int nbNonZero){

	for (int i=0; i<nbNonZero; i++)
	{
		for (int j=0; j<3; j++)
			printf("%f ", sparceA[3*i+j]);

		printf("\n");
	}
}
