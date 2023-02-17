#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include <omp.h>
#include <libgen.h>
#include "utils/hashtable.h"
#include "utils/blas.h"
#include "utils/utils.h"
#include "page_rank.h"

// ******* MAIN *******
int main(int argc, char **argv) {
    double temps, debut, fin;
    unsigned char *adjacency_matrix;
    int matrix_dim, 
        nbNonZero = -1,
        *sparce_adjacency_matrix,
        i;

    if(argc != 5){
        fprintf(stderr, "Usage: %s path_adjancecy_matrix_file parallel<0|1> sparce<0|1> outputFiles<0|1>\n", argv[0]);
        exit(1);
    }

    // Parse arguments
    char * path_adjancecy_matrix_file = argv[1];
    int MODE_EXEC = atoi(argv[2]);
    const int SPARCE_REPRESENTATION = atoi(argv[3]);
    int outputfiles = atoi(argv[4]);

    // Get directory of the adjancecy file
    char* directory_path = dirname(strdup(path_adjancecy_matrix_file));
    printf("directory path : %s\n", directory_path);

    // *************** BEGIN : READ adjacency matrix **********************
    debut = omp_get_wtime();

    nbNonZero = read_adj_matrix_from_file(path_adjancecy_matrix_file, &adjacency_matrix, &sparce_adjacency_matrix, &matrix_dim, SPARCE_REPRESENTATION);
    
    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create adjacency_matrix ==>  : %lf s\n",temps);
    // *************** END : CREATE adjacency matrix from file **********************

    // ************************ BEGIN : CREATE TRANSITION MATRIX **************************************
    debut = omp_get_wtime();

    double* transition_matrix = create_transition_matrix(adjacency_matrix, sparce_adjacency_matrix, matrix_dim, nbNonZero, SPARCE_REPRESENTATION, MODE_EXEC);
    
    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create transition_matrix ==>  : %lf s\n",temps);
    // ************************ END : CREATE TRANSITION MATRIX **************************************

    // ****** BEGIN : CREATE RANDOM TELEPORTATION VECTOR ******
    double *teleportation_vector = malloc(matrix_dim*sizeof(double));
    {
        #pragma omp parallel for schedule(static) private(i)
        for (i = 0; i < matrix_dim; i++)
        {
            teleportation_vector[i] = 1./matrix_dim;
        }    
    }
    // ****** END : CREATE RANDOM TELEPORTATION VECTOR ******


    // ************************************** BEGIN :  Apply PageRank Algorithm **************************************
	debut = omp_get_wtime();
    double *pg_vector = page_rank(transition_matrix, teleportation_vector, nbNonZero, matrix_dim, 0.85, 0.00001, MODE_EXEC, SPARCE_REPRESENTATION);
	fin = omp_get_wtime();
	temps = fin - debut;

	printf("Time PageRank Algorithm  ==>  : %lf s\n",temps);
    // ************************************** END :  Apply PageRank Algorithm **************************************


    // ******* BEGIN : Write Results to DISK ********
    // WRITE TRANSITION MATRIX TO DISK
    if(outputfiles){
        char transition_matrix_file_path[254] = "";
        strcpy(transition_matrix_file_path, directory_path);
        strcat(transition_matrix_file_path, "/transition_matrix.csv");

        write_matrix(matrix_dim, matrix_dim, nbNonZero, transition_matrix, transition_matrix_file_path, SPARCE_REPRESENTATION);
    }

    // WRITE PAGERANK VECTOR TO DISK
    char pagerank_vector_file_path[254] = "";
    strcpy(pagerank_vector_file_path, directory_path);
    strcat(pagerank_vector_file_path, "/pagerank_vector.csv");
    write_matrix(1, matrix_dim, nbNonZero, pg_vector, pagerank_vector_file_path, 0);
    // ******* END : Write Results to DISK ********

    // Free allocated space
    if(SPARCE_REPRESENTATION)
        free(sparce_adjacency_matrix);
    else
        free(adjacency_matrix);


    free(transition_matrix);
    free(pg_vector);

    return 0;
}

