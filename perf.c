#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "hashtable.h"
#include "blas.h"
#include "page_rank.h"
#include <omp.h>
#include <libgen.h>

// ******* BEGIN : Declarations *******

void write_adjacency_matrix(int row, int col, double* matrix, char* filepath);
void write_matrix(int row, int col, double* matrix, char* filepath);
double* create_transition_matrix(double* adjacency_matrix, double *out_links_vector, int n, int nbNonZero, int sparce_rep, int parralel);
double* read_adj_matrix_from_file(int nb_ligne, int nb_colone, double* matrix, FILE *matrix_file, int* nb_non_zero);

// ******* END  : Declarations *******

// ******* MAIN *******
int main(int argc, char **argv) {
    double temps, debut,fin,
            *opposite_rep_adjacency_matrix,
            *adjacency_matrix;
    int matrix_dim, nbNonZero = -1;

    if(argc < 6){
        fprintf(stderr, "Usage: %s path_adjancecy_matrix_file matrix_dim parallel<0|1> sparce<0|1> outputFiles<0|1>\n", argv[0]);
        exit(1);
    }

    char * path_adjancecy_matrix_file = argv[1];
    printf("path_adjancecy_matrix_file : %s\n", path_adjancecy_matrix_file);


    // Get directory of the adjancecy file
    char* directory_path = dirname(strdup(path_adjancecy_matrix_file));
    printf("directory path : %s\n", directory_path);

    matrix_dim = atoi(argv[2]);
    int MODE_EXEC = atoi(argv[3]);
    int SPARCE_REPRESENTATION = atoi(argv[4]);
    int outputfiles = atoi(argv[5]);

    // SO(n^2)
    // Load adjacency matrix from file
    debut = omp_get_wtime();

    adjacency_matrix = calloc(matrix_dim * matrix_dim, sizeof(double));
    FILE* mat = fopen(path_adjancecy_matrix_file, "r");
    if(mat == NULL) {
        perror("Error opening file");
        exit(1);
    }
    read_adj_matrix_from_file(matrix_dim, matrix_dim, adjacency_matrix, mat, &nbNonZero);
    fclose(mat);

    // Create sparce matrix representation of the adjacency matrix
    if(SPARCE_REPRESENTATION){
        opposite_rep_adjacency_matrix = adjacency_matrix;
        // TO(n^2), SO(nzero * 3)
        adjacency_matrix = matrix_to_sparce(adjacency_matrix, matrix_dim, matrix_dim, &nbNonZero);
    }


    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create adjacency_matrix ==>  : %lf s\n",temps);


    // ************************ BEGIN : CREATE TRANSITION MATRIX **************************************
    debut = omp_get_wtime();
    // SO(n)
    double *out_links_vector = calloc(matrix_dim,sizeof(double));
    // SO(n)
    double *vector_of_ones = calloc(matrix_dim,sizeof(double));

    // TO(n / p)
    #pragma omp parallel for schedule(static) if(MODE_EXEC)
    int i;
for(i = 0; i < matrix_dim; i++){
        vector_of_ones[i] = 1;
    }

    // Get out_links_vectors
    if(SPARCE_REPRESENTATION)
        /*
            Time complexity : TO(nzero)
        */
        Sparce_Matrix_Vector_Product(adjacency_matrix, vector_of_ones, nbNonZero, matrix_dim, out_links_vector, MODE_EXEC);
    else
        /*
            Sequential case: TO(n^2), SO(1)
            Parallel case: TO(n^2 / p), SO(1)
        */
        Matrix_Vector_Product(adjacency_matrix, vector_of_ones, matrix_dim, matrix_dim, out_links_vector, MODE_EXEC);
        
    free(vector_of_ones);

    debut = omp_get_wtime();
    /*
        Time/Space complexity :
        Sequential case:
            using normal representation: TO(n^2) SO(n^2)
            using sparce representation: TO(nzero) SO(nzero * 3)
        Parallel case:
            using normal representation: TO(n^2 / p) SO(n^2)
            using sparce representation: TO(nzero / p) SO(nzero * 3)
    */
    double* transition_matrix = create_transition_matrix(adjacency_matrix, out_links_vector, matrix_dim, nbNonZero, SPARCE_REPRESENTATION, MODE_EXEC);
    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create transition_matrix ==>  : %lf s\n",temps);
    // ************************ END : CREATE TRANSITION MATRIX **************************************


    // Apply PageRank Algorithm
    /*
        Time/Space complexity :
        Sequential case:
            Normal representation: TO(nb_iterations * n^2) SO(4n)
            Sparce representation: TO(nb_iterations * nzero) SO(3n)
        Parallel case:
            Normal representation: TO(nb_iterations * (n^2 / p)) SO(4n)
            Sparce representation: TO(nb_iterations * (nzero + n/p)) SO(4n)
    */
	debut = omp_get_wtime();
    double *pg_vector = page_rank(transition_matrix, nbNonZero, matrix_dim, 0.85, 0.00001, MODE_EXEC, SPARCE_REPRESENTATION);
	fin = omp_get_wtime();
	temps = fin - debut;

	printf("Time PageRank Algorithm  ==>  : %lf s\n",temps);

    // WRITE TRANSITION MATRIX TO DISK
    if(outputfiles){
        char transition_matrix_file_path[254] = "";
        strcpy(transition_matrix_file_path, directory_path);
        strcat(transition_matrix_file_path, "/transition_matrix.csv");

        // @TODO free sparce_to_matrix
        if(SPARCE_REPRESENTATION)
            write_matrix(matrix_dim, matrix_dim, sparce_to_matrix(transition_matrix, NULL, matrix_dim, matrix_dim, nbNonZero), transition_matrix_file_path);
        else
            write_matrix(matrix_dim, matrix_dim, transition_matrix, transition_matrix_file_path);
    }

    // WRITE PAGERANK VECTOR TO DISK
    char pagerank_vector_file_path[254] = "";
    strcpy(pagerank_vector_file_path, directory_path);
    strcat(pagerank_vector_file_path, "/pagerank_vector.csv");
    write_matrix(1, matrix_dim, pg_vector, pagerank_vector_file_path);


    // Free allocated space
    free(out_links_vector);
    free(adjacency_matrix);
    free(transition_matrix);

    return 0;
}


double* read_adj_matrix_from_file(int nb_ligne, int nb_colone, double* matrix, FILE *matrix_file, int* nb_non_zero)
{
    *nb_non_zero = 0;
    int i;
for(i = 0; i < nb_ligne; i++)
    {
        int j;
for(j = 0; j < nb_colone; j++)
        {
            double v = 0;
            fscanf(matrix_file, "%lf", &v);
            matrix[i * nb_ligne + j] = v;
            if(v != 0)
                *nb_non_zero = *nb_non_zero + 1;
        }
    }
    return matrix;
}

// Write  matrix to a file
void write_matrix(int row, int col, double* matrix, char* filepath){
    FILE* stream = fopen(filepath, "w");
    if(stream == NULL) {
        perror("Error opening file");
        exit(1);
    }
    int i;
for(i = 0; i < row; i++){
        int j;
for(j = 0; j < col; j++){
            fprintf(stream,"%0.4f ", matrix[i*col+j]);
        }
        fprintf(stream,"\n");
    }

    fclose(stream);
}

/*
    n*n : size adjacency matrix
    nzero : size sparce matrix (i.e number of non-zero elements of adjacency matrix)

    Time/Space complexity :
    Sequential case:
        using normal representation: TO(n^2) SO(n^2)
        using sparce representation: TO(nzero) SO(nzero * 3)
    Parallel case:
        using normal representation: TO(n^2 / p) SO(n^2)
        using sparce representation: TO(nzero / p) SO(nzero * 3)
*/
double* create_transition_matrix(double* adjacency_matrix, double *out_links_vector, int n, int nbNonZero, int sparce_rep, int parralel){
    // Initialize transition matrix
    double *transition_matrix;
    if(sparce_rep){
        // SO(nbNonZero * 3)
        transition_matrix = calloc(nbNonZero * 3, sizeof(double));

        // TO(nzero / p)
        #pragma omp parallel for schedule(static) if(parralel)
        int i;
for(i = 0; i < nbNonZero; i++){
                // if link exists from j to i
                // then probabilty of going to i from j is 1 div number of out links of node j
                int j = adjacency_matrix[i*3];
                transition_matrix[i*3] = adjacency_matrix[i*3+1];
                transition_matrix[i*3+1] = j;
                transition_matrix[i*3+2] = 1 / out_links_vector[j];
        }
    }
    else{
        // SO(n^2)
        transition_matrix = calloc(n * n, sizeof(double));

        // TO(n^2 / p)
        #pragma omp parallel for schedule(static) if(parralel)
        int i;
for(i = 0; i < n; i++){
            int j;
for(j = 0; j < n; j++){
                // if link exists from j to i
                // then probabilty of going to i from j is 1 div number of out links of node j
                if(adjacency_matrix[j*n+i])
                    transition_matrix[i*n+j] = 1 / out_links_vector[j];
            }
        }
    }

    return transition_matrix;
}
