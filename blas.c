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
    int i;
    for (i = 0; i < n; i++){
        x[i] = (rand() % (min - max + 1)) + min;
    }
    return x;
}

double* initRandomMatrix(double* A, int row,int col, int min, int max){
    if(A == NULL)
        A = calloc(row*col, sizeof(double));
    int i;
for(i = 0; i < row; i++){
        int j;
for(j = 0; j < col; j++){
            A[i*col+j] = (rand() % (min - max + 1)) + min;
        }
    }
    return A;
}

void displayVector(double* x, int n){
    printf("\n");
    int i;
for(i = 0; i < n; i++){
            printf("%f ", x[i]);
    }
    printf("\n");
}

void displayMatrix(int row, int col, double* matrix){
    printf("_______________________________\n");
    int i;
for(i = 0; i < row; i++){
        int j;
for(j = 0; j < col; j++){
            printf("%0.2f  ", matrix[i*col+j]);
        }
        printf("\n");
    }
    printf("_______________________________\n");
}

double* transposeMatrix(int row, int col, double* A){
    double* AT = calloc(row*col, sizeof(double));
    int i;
for(i=0; i<row; i++) {
        int j;
for(j=0; j<col; j++) {
            AT[j*row+i] = A[i*col+j];
        }
    }
    return AT;
}

void Vector_Vector_Addition(double* x, double* y, int n){
    int i;
for(i = 0; i < n; i++){
        y[i] = x[i] + y[i];
    }
}

// Produit scalaire vecteur
void Vector_Scalar_Product(double* x, double alpha, int n) {
    int i;
for(i = 0; i < n; i++){
        x[i] *=  alpha;
    }
}

// Produit scalaire vecteur
void Vector_Scalar_Product_parallel(double* x, double alpha, int n) {
    int i;
    #pragma omp parallel for schedule(static) private(i)
    for(i = 0; i < n; i++){
        x[i] *=  alpha;
    }
}


/* 
    Produit scalaire 2 vecteurs
    p : nb processors
    Time/Space complexity :
        Sequential case: TO(n) SO(1)
        Parallel case: TO(n / p) SO(1)
*/
double DotProduct(double* x, double* y, int n, int parallel) {
    double result = 0;
    int i;
    #pragma omp parallel for schedule(static) reduction(+:result) private(i) if(parallel)
    for(i = 0; i < n; i++){
        result += x[i] * y[i];
    }
    return result;
}

/*
    n*n : size matrix A
    Time/Space complexity: 
        Sequential case: TO(n^2), SO(1)
        Parallel case: TO(n^2 / p), SO(1)
*/
void Matrix_Vector_Product(double* A, double* v, int row, int col, double* Av, int parallel){
    if(parallel)
        Matrix_Vector_Product_parralel(A, v, row, col, Av);
    else
        Matrix_Vector_Product_sequential(A, v, row, col, Av);
}

void Matrix_Vector_Product_uchar(unsigned char* A, double* v, int row, int col, double* Av, int parallel){
    if(parallel)
        Matrix_Vector_Product_parralel_uchar(A, v, row, col, Av);
    else
        Matrix_Vector_Product_sequential_uchar(A, v, row, col, Av);
}

/*
    n*n : size matrix A
    Time/Space complexity: TO(n^2), SO(1)
*/
void Matrix_Vector_Product_sequential(double* A, double* v, int row, int col, double* Av){
    int i;
for(i=0; i<row; i++) {
        Av[i] = 0;
        int j;
for(j=0; j<col; j++) {
            Av[i] += A[i*col+j] * v[j];
        }
    }
}

void Matrix_Vector_Product_sequential_uchar(unsigned char* A, double* v, int row, int col, double* Av){
    int i;
    for(i=0; i<row; i++) {
        Av[i] = 0;
        int j;
        for(j=0; j<col; j++) {
            Av[i] += A[i*col+j] * v[j];
        }
    }
}

/*
    n*n : size matrix A
    Time/Space complexity: TO(n^2 / p), SO(1)
*/
void Matrix_Vector_Product_parralel(double* A, double* v, int row, int col, double* Av){
    omp_set_nested(2);
    int i;
    #pragma omp parallel for schedule(static) private(i)
    for(i = 0; i<row; i++){
        double sum =  0;
        int j;
        for(j=0; j<col; j++) {
            sum += A[i*col+j] * v[j];
        }
        Av[i] = sum;
    }
}

void Matrix_Vector_Product_parralel_uchar(unsigned char* A, double* v, int row, int col, double* Av){
    omp_set_nested(2);
    int i;
    #pragma omp parallel for schedule(static) private(i)
    for(i = 0; i<row; i++){
        double sum =  0;
        int j;
        for(j=0; j<col; j++) {
            sum += A[i*col+j] * v[j];
        }
        Av[i] = sum;
    }
}



/*
    vres = α.Ax+βy
    n*n : size matrix A
    nzero : size sparce matrix (i.e number of non-zero elements of matrix A)

    Time/Space complexity :
    Sequential case:
        using normal representation: TO(n^2) SO(n)
        using sparce representation: TO(nzero + n) SO(n)
    Parallel case:
        using normal representation: TO(n^2 / p) SO(n)
        using sparce representation: TO(nzero / p + n + n/p) SO(n*p + n)
*/
void blas21(double* A, double* x, double* y, double* vres, double alpha, double beta, int row, int col, int nbNonZeroA, int parallel, int sparce_rep){
    if(parallel){
        if(!sparce_rep)
            blas21_parallel(A, x, y, vres, alpha, beta, row, col);
        else
            blas21_parallel_sparce(A, x, y, vres, alpha, beta, col, nbNonZeroA);
    }
    else{
        if(!sparce_rep)
            blas21_sequential(A, x, y, vres, alpha, beta, row, col);
        else
            blas21_sequential_sparce(A, x, y, vres, alpha, beta, col, nbNonZeroA);
    }
}


/*
    x= α.Ax+βy
    n : matrix size
    Time complexity : TO(n^2)
    Space complexity : SO(n)
*/
void blas21_sequential(double* A, double* x, double* y, double* vres, double alpha, double beta, int row, int col){
    // double* v = malloc(row * sizeof(double));
    // TO(2n^2 + 3n)
    int i;
for(i=0; i<row; i++) {
        double ax = 0;
        int j;
for(j=0; j<col; j++) {
            ax += A[i*col+j] * x[j];
        }
        vres[i] = alpha * ax + beta * y[i];
    }
    // TO(n)
    // memcpy(x, v, row * sizeof(double));    
    // free(v);
}

/*
    x= α.Ax+βy
    n*n : size matrix
    nzero : size sparce matrix (i.e number of non-zero elements)
    Time complexity : TO(nzero + n)
    Space complexity : SO(n)
*/
void blas21_sequential_sparce(double* sparceA, double* x, double* y, double* vres, double alpha, double beta, int n, int sizeSparceA){
    // SO(n)
    double* Av = calloc(n ,sizeof(double));
    // Matrix Vector Product 
    Sparce_Matrix_Vector_Product(sparceA, x, sizeSparceA, n, Av, 0);
    // TO(n)
    int j;
for(j=0; j<n; j++) {
        vres[j] = alpha * Av[j] + beta * y[j];
    }
    free(Av);
}

/*
    vres= α.Ax+βy
    n : matrix size
    Time complexity : TO(n^2 / p)
    Space complexity : SO(n)
*/
void blas21_parallel(double* A, double* x, double* y, double* vres, double alpha, double beta, int row, int col){
    // SO(n)
    // double* v = malloc(row * sizeof(double));
    int i;
    #pragma omp parallel for schedule(static) private(i)
    // TO(n^2 / p)
    for(i = 0; i<row; i++) {
        double ax = 0;
        int j;
for(j=0; j<col; j++) {
            ax += A[i*col+j] * x[j];
        }
        vres[i] = alpha * ax + beta * y[i];
    }

    // memcpy(x, v, row * sizeof(double));    
    // free(v);
}

/*
    vres = α.Ax+βy
    p : nb processors
    n*n : size matrix
    nzero : size sparce matrix (i.e number of non-zero elements)
    Time complexity : TO(nzero + n/p)
    Space complexity : SO(n)
*/
void blas21_parallel_sparce(double* sparceA, double* x, double* y, double* vres, double alpha, double beta, int n, int sizeSparceA){
    // SO(n)
    double* Av = calloc(n ,sizeof(double));
    /*
        Matrix Vector Product : TO(nzero)
    */
    Sparce_Matrix_Vector_Product(sparceA, x, sizeSparceA, n, Av, 1);
    // TO(n / p)
    int j;
    #pragma omp parallel for schedule(static) private(j)
    for( j=0; j<n; j++) {
        vres[j] = alpha * Av[j] + beta * y[j];
    }
    free(Av);
}


double* Matrix_Matrix_Product(double* A, double* B, int rowA, int colA, int colB, double* AB){
    if(AB == NULL)
        AB = calloc(rowA*colB, sizeof(double));
    int i;
for(i=0; i<rowA; i++) {
        int k;
for(k=0; k<colB; k++) {
            int j;
for(j=0; j<colA; j++) {
                AB[i*colB+k] += A[i*colA+j] * B[j*colB+k];
            }
        }
    }
    return AB;
}


double* Matrix_Matrix_Subsctraction(double* A, double* B, int row, int col, double* AB){
    if(AB == NULL)
        AB = calloc(row*col, sizeof(double));
    int i;
for(i=0; i<row; i++) {
        int k;
for(k=0; k<col; k++) {
            AB[i*col+k] = A[i*col+k] - B[i*col+k];
        }
    }
    return AB;
}

/*
    Norme2 d'un vector
    Time complexity : O(2n + 1)
    Space complexity : O(1)
*/
double Norme(double* x, int n, int parallel){
    // Time complexity : O(2n + 1)
    // Space complexity : O(1)
    return sqrt(DotProduct(x, x, n, parallel));
}

// Norme1 d'un vector
double Norme_One(double* x, int n){
    double result = 0;
    int i;
for(i = 0; i < n; i++){
        result += x[i];
    }
    return result;
}

// Norme1 d'un vector
double Norme_One_parallel(double* x, int n){
    double result = 0;
    int i;
    #pragma omp parallel for schedule(static) reduction(+:result) private(i)
    for(i = 0; i < n; i++){
        result += x[i];
    }
    return result;
}

// La norme Frobenius d’une matrice
double NormeFrobenius(int nb_ligne, int nb_colonne, double* matrice){
    double result = 0;
    int i;
for(i=0; i<nb_ligne; i++) {
        int j;
for(j=0; j<nb_colonne; j++) {
            result += pow(matrice[i*nb_colonne+j], 2);
        }
    }
    return sqrt(result);
} 



/*   
    n*n : size square matrix
    nzero : size sparce matrix (i.e number of non-zero elements in the original matrix)
    Sequential case:
        Time complexity :  TO(nzero)
*/
void Sparce_Matrix_Vector_Product(double* A, double* v, int sizeA, int sizeV, double* Av, int parallel){

        int i;
        for(i=0; i<sizeA; i++) {
            int r_idx = A[i*3];
            int c_idx = A[i*3+1];
            Av[r_idx] += A[i*3+2] * v[c_idx];
        }
}

void Sparce_Matrix_Vector_Product_int(int* A, double* v, int sizeA, int sizeV, double* Av, int parallel){

        int i;
        for(i=0; i<sizeA; i++) {
            int r_idx = A[i*3];
            int c_idx = A[i*3+1];
            Av[r_idx] += A[i*3+2] * v[c_idx];
        }
}



/*
    Count number of non-zero elements of matrix A
    n*n : size square matrix A
    p : number of processors
    Time/ Space complexity:
        Sequential case: TO(n^2), SO(1)
        Parallel case: TO(n^2 / p)
*/
int count_zero_matrix(double *A, int row, int col){
	int nbNonZero = 0;
    int i;
for(i=0; i<row; i++) {
        int j;
for(j=0; j<col; j++) {
            if (A[i*col+j] != 0)
				nbNonZero++;
        }
    }
    return nbNonZero;
}

/* 
    Convert normal matrix representation to sparce matrix representation
    Time/Space complexity:
        sequential case: TO(n^2), SO(nzero * 3)
        parallel case:
*/
double* matrix_to_sparce(double *A, int row, int col, int* nbNonZero){
	// If nbNonZero not known
	if(*nbNonZero < 0){
		*nbNonZero = count_zero_matrix(A, row, col);
	}

	double *sparceA = calloc(*nbNonZero * 3,sizeof(double)); // SO(nzero * 3)

	// Making of new matrix 
	int k = 0;
    // TO(n^2)
	int i,j;
    for(i = 0; i < row; i++)
        for(j = 0; j < col; j++)
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
	int i;
for(i=0; i<nbNonZero; i++)
	{
		r_idx = sparceA[i*3];
		c_idx = sparceA[i*3+1];
		A[r_idx*col+c_idx] = sparceA[3*i+2];
	}

	return A;
}

void display_sparce_matrix(double *sparceA, int nbNonZero){

	int i;
for(i=0; i<nbNonZero; i++)
	{
		int j;
for(j=0; j<3; j++)
			printf("%f ", sparceA[3*i+j]);

		printf("\n");
	}
}
