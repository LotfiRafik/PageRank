/*
    Page Rank Algorithm 

    n*n : size square matrix 
    nzero : number of non-zero elements
    Time/Space complexity :
        Sequential case:
            Normal representation: TO(nb_iterations * n^2) SO(4n)
            Sparce representation: TO(nb_iterations * nzero) SO(4n)
        Parallel case:
            Normal representation: TO(nb_iterations * (n^2 / p)) SO(4n)
            Sparce representation: TO(nb_iterations * n) SO(n*p + 4n)
*/
double *page_rank(double *A, int nbNonZeroA , int n, double B, double p, int parallel_mode, int sparce_rep);