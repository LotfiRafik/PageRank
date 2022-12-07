#include <stdio.h>
#include <stdlib.h>
#include "blas.h"
#include "page_rank.h"
#include <string.h>

double *page_rank(double *A, int n, double B, double p, int parallel_mode){

    double *x=calloc(n , sizeof(double)),
            *x_prec=calloc(n , sizeof(double)),
             *teleportation_vector=calloc(n , sizeof(double)),
             error;

    /* TODO 
        create parallel function to initialize a vector with a scalar value
    */
    // initializer le vecteur de depart
    for (int j=0; j<n; j++){
        x[j]=1./n;
        // TODO Personalized pageRank changes this vector
        teleportation_vector[j] = 1./n;
    }

    // i represente le nombre d'itirations
    int  nb_iterations = 0;
    
    // boucler jusqu'à la convergence
    do{
        // incrémenter le nombre d'itiration
        nb_iterations++;

        // sauvgarder le vecteur i-1
        memcpy(x_prec, x, n * sizeof(double));    

        // TODO Parallelize the for loop
        blas21(A, x, teleportation_vector, B, 1-B, n, n, parallel_mode);

        /* TODO 
            Check with quentin if we use norme 1 or norme 2
            in the excercice TD quentin did not normalize the vector x
        */
        // normaliser le vecteur x
        double norm = Norme(x, n);
        // Vector_Scalar_Product(x, 1/norm, n);

        // calculer l'erreur entre le nouveau vecteur et le vecteur précédent
        double y[n]; 
        // TODO Parralelize the for loop
        for (int i = 0 ; i<n; i++){
           y[i] = x[i] - x_prec[i]; 
        }
        // verifier la condition de convergence
        error = Norme(y, n);

    } while(error > p);


    printf("nombre d'iteration: %d \n", nb_iterations);

    // Free allocated memory space
    free(x_prec);
    free(teleportation_vector);
    return x;
}