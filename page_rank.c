#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "utils/blas.h"
#include "page_rank.h"


/*
    Page Rank Algorithm 

    n*n : size square matrix 
    nzero : number of non-zero elements
    Sequential case:
        Using normal representation:
            Time complexity : O(nb_iterations * n^2)
            Space complexity : O(4n)
        Using sparce representation:
            Time complexity : O(nb_iterations * nzero)
            Space complexity : O(4n)
    Parallel case:
        Using normal representation:
            Time complexity : O(nb_iterations * (n^2 / p))
            Space complexity : O(4n)
        Using sparce representation:
            Time complexity : TO(nb_iterations * (nzero + n/p)) 
            Space complexity : SO(4n)
*/
double *page_rank(double *A, double* teleportation_vector, int nbNonZeroA , int n, double B, double p, int parallel_mode, int sparce_rep){

    // Space complexity : Auxiliary + Space = O(3n + 1)
    double *x=calloc(n , sizeof(double)),
            *x_prec=calloc(n , sizeof(double)),
            error;

    int sum_teleportation_vector = 0;
    int i;

    srand(time(NULL));   // Initialization, should only be called once.

    #pragma omp parallel for schedule(static) private(i)
    for (i = 0; i < n; i++)
    {
        // initializer le vecteur de depart
        x[i] = 1./n;
    }


    // nombre d'iterations
    int  nb_iterations = 0;
    
    // boucler jusqu'à la convergence
    do{
        // incrémenter le nombre d'itiration
        nb_iterations++;

        // sauvgarder le vecteur i-1
        memcpy(x_prec, x, n * sizeof(double));    // TO(n)

        /*
            x= α.Ax+βy
            n*n : size matrix A
            nzero : size sparce matrix (i.e number of non-zero elements of matrix A)

            Time/Space complexity :
            Sequential case:
                using normal representation: TO(n^2) SO(n)
                using sparce representation: TO(nzero + n) SO(n)
            Parallel case:
                using normal representation: TO(n^2 / p) SO(n)
                using sparce representation: TO(nzero + n/p) SO(n*p + n)
        */
        blas21(A, x_prec, teleportation_vector, x, B, 1-B, n, n, nbNonZeroA, parallel_mode, sparce_rep);

        /* TODO 
            Check with quentin if we use norme 1 or norme 2
            in the excercice TD quentin did not normalize the vector x
        */
        // normaliser le vecteur x
        double norm = Norme(x, n, parallel_mode);
        // Vector_Scalar_Product(x, 1/norm, n);

        /* 
            calculer l'erreur entre le nouveau vecteur et le vecteur précédent
            Space complexity : Auxiliary + Space = O(n)
        */
        double y[n]; 
        // TO(n / p)
        #pragma omp parallel for schedule(static) private(i) if(parallel_mode)
        for(i = 0; i<n; i++){
           y[i] = x[i] - x_prec[i]; 
        }

        // verifier la condition de convergence
        error = 0;
        // calculer la norme : ||x[i] - x_prec[i]||
        #pragma omp parallel for schedule(static) reduction(+:error) private(i) if(parallel_mode)
        for(i = 0; i < n; i++){
            double y = x[i] - x_prec[i];
            error += y * y;
        }
        error = sqrt(error);


        // error = Norme(y, n, parallel_mode);

    } while(error > p);


    printf("nombre d'iteration: %d \n", nb_iterations);

    // Free allocated memory space
    free(x_prec);
    return x;
}