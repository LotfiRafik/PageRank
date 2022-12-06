#include <stdio.h>
#include <stdlib.h>
#include "blas.h"


/*
=================================================================
Lecture de la matrice

*/
double *Lecture_Matrice(int nb_ligne, int nb_colone, FILE *matrice)
{

    double *mat = malloc(nb_ligne * nb_colone * sizeof(double));

    for (int i = 0; i < nb_ligne; i++)
    {
        for (int j = 0; j < nb_colone; j++)
        {
            double v = 0;
            fscanf(matrice, "%lf", &v);

            mat[j * nb_ligne + i] = v;
        }
    }

    return mat;
}

/*
===========================================================================
affichage de la matrice

*/




double *page_rank(double *A, int n, double B, double p){

    double *x=calloc(n , sizeof(double));

    double *x_new=calloc(n , sizeof(double));


    // initializer le vecteur de depart
    for (int j=0; j<n; j++){
        x_new[j]=1./n;
    }


    // transposer la matrice pour le calcul de vecteur matrice
    double *Z = transposeMatrix(n,n,A);


    // i represente le nombre d'itirations
    int  i = 0;
    
  

    
    // boucler jusqu'à la convergence
    for (;;){


        // incrémenter le nombre d'itiration
        i++;

        // produit matrice vecteur
        Matrix_Vector_Product(Z, x_new, n,n , x);


        // claculer le vecteur Page Rank avec le coefficient d'amortissement
        for(int i = 0; i<n; i++){
            x[i]=B * x[i]+(1-B)/n;
        }

        // normaliser le vecteur x
        double norm = Norme(x, n);
        for(int i = 0; i<n; i++){
            x[i]= x[i]/norm;
        }


        // claculer l'erreur entre le nouveau vecteur et le vecteur précédent
        double y[n]; 

        for ( int i = 0 ; i<n; i++){
           y[i] = x[i] - x_new[i]; 
        }

        // verifier la condition de convergence
        double error = Norme(y, n );
        if ( error < p ){
            break;
        }

        // initializer le vecteur initial
        for ( int i = 0 ; i<n; i++){
            x_new[i] = x[i]; 
        }
   


        
    }

    printf("nombre d'iteration: %d \n", i);

    return x;


    free(x_new);
    free(Z);



    
}



int main(int argc, char* argv){
    FILE *mat;
    mat = fopen("Deezer-small-DS/Exo Td/transition_matrix.csv", "r");
    if (mat == NULL)
    {
        printf("le fichier n'existe pas");
        return 1;
    }

    double *mat1 = Lecture_Matrice(4, 4, mat);
    fclose(mat);



    double *x=page_rank(mat1, 4, 0.85, 0.001);

    displayVector(x, 4);

    free(mat1);
    free(x);

    

   
    return 0;
}