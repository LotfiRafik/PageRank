#include <stdio.h>
#include <stdlib.h>
#include <blas.h>


double *page_rank(double *A, int n, double B, double p){

    double *x=malloc(n * sizeof(double));

    double *x_new=malloc(n * sizeof(double));

    for (int j=0; j<n; j++){
        x_new[j]=1./n;
    }

    int i =0;

    for (;;){
        
        Matrix_Vector_Product(A, x_new, n,n , x);


        for(int i = 0; i<n; i++){
            x[i]=B * x[i] + (1-B)/n ;
        }
        double norm = Norme(x, n);

        for(int i = 0; i<n; i++){
            x[i]= x[i]/norm;
        }

        double y[n]; 

        for ( int i = 0 ; i<n; i++){
           y[i] = x[i] - x_new[i]; 
        }


        double error = Norme(y, n );


        if ( error < p ){
            return x;
        }

        for ( int i = 0 ; i<n; i++){
            x_new[i] = x[i]; 
        }


        
    }

   

    return x;


    
}



int main(int argc, char* argv){
   
    return 0;
}