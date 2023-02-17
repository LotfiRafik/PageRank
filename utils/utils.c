#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include "utils.h"
#include "blas.h"

// Separates the string into substrings, splitting the string into substrings 
// based on the separator characters (i.e separators).  The function returns an
// array of pointers to strings, dynamically allocated on the heap, and it 
// effectively "returns" the number of these strings via pass-by-pointer using 
// the parameter count.  
// Time complexity : O(2n) with n being length of string to split
char **split(char *string, char seperator, int* nb_tokens){
    // get the length of the string
    int len = strlen(string);
    int i = 0 ,
        j = 0;

    *nb_tokens = 0;

    // First loop to find out how many tokens there is
    while (i <= len){
        j = 0;
        while(i < len && string[i] != seperator){
            i++;
            j++;
        }
        // ****** Token found *****
        *nb_tokens = *nb_tokens + 1;

        // skip current separator
        i++;
    }

    
    // Maximum token's length is the len of the string (i.e string contains no separators)
    char token[len];
    
    // allocate space for a dynamically allocated array of *nb_tokens* number of 
    // pointers to strings
    char **tokens = malloc(sizeof(char*) * *nb_tokens);
    
    // Reinitialize count variables
    *nb_tokens = 0;
    i = 0;
    j = 0;

    // Second loop to extract tokens
    while (i <= len){
        // Initialize token
        token[0] = '\0';
        
        j = 0;
        while(i < len && string[i] != seperator){
            token[j] = string[i];
            i++;
            j++;
        }
        // ****** Token found *****

        // add a null terminator on to the end of token to terminate the string
        token[j] = '\0';

        tokens[*nb_tokens] = strdup(token);

        *nb_tokens = *nb_tokens + 1;

        // skip current separator
        i++;
    }
    
    // return our array of strings  
    return tokens;
}

void free_splitted_tokens(char** tokens, int size){
    int i;
for(i=0; i<size; i++)
        free(tokens[i]);
        
    free(tokens);
}

// Write  matrix to a file
void write_adjacency_matrix(int row, int col, int nbNonZero, unsigned char* adjacency_matrix, int* sparce_adjacency_matrix, char* filepath, int sparce_rep){
    int i,j;

    FILE* stream = fopen(filepath, "w");
    if(stream == NULL) {
        perror("Error opening file");
        exit(1);
    }
    // If matrix is stored in compressed representation
    if(sparce_rep){
        fprintf(stream, "%d %d %d\n", row, col, nbNonZero);
        for (i = 0; i < nbNonZero; i++)
        {
            fprintf(stream,"%d %d %d\n", sparce_adjacency_matrix[i*3], sparce_adjacency_matrix[i*3+1], 1);
        }
    }
    else{
        fprintf(stream, "%d %d %d\n", row, col, row * col);
        for(i = 0; i < row; i++)
            for(j = 0; j < col; j++)
                fprintf(stream,"%d %d %u\n", i, j, adjacency_matrix[i*col+j]);
    }

    fclose(stream);
}

// Write  matrix to a file
void write_matrix(int row, int col, int nbNonZero, double* matrix, char* filepath, int sparce_rep){
    int i,j;

    FILE* stream = fopen(filepath, "w");
    if(stream == NULL) {
        perror("Error opening file");
        exit(1);
    }
    // If matrix is stored in compressed representation
    if(sparce_rep){
        fprintf(stream, "%d %d %d\n", row, col, nbNonZero);
        for (i = 0; i < nbNonZero; i++)
        {
            fprintf(stream,"%d %d %0.4lf\n", (int)matrix[i*3], (int)matrix[i*3+1], matrix[i*3+2]);
        }
    }
    else{
        fprintf(stream, "%d %d %d\n", row, col, row * col);
        for(i = 0; i < row; i++)
            for(j = 0; j < col; j++)
                fprintf(stream,"%d %d %0.4lf\n", i, j, matrix[i*col+j]);
    }

    fclose(stream);
}

// Write vector to a file
void write_vector(int dim, double* vector, char* filepath){
    int i,j;

    FILE* stream = fopen(filepath, "w");
    if(stream == NULL) {
        perror("Error opening file");
        exit(1);
    }
    
    for(i = 0; i < dim; i++)
        fprintf(stream,"%0.4lf\n", vector[i]);

    fclose(stream);
}


int read_adj_matrix_from_file(char* path_adjancecy_matrix_file, unsigned char** adjacency_matrix, int** sparce_adjacency_matrix, int *matrix_dim, int sparce_rep)
{
    char line[1024];
    int nb_line, nb_col, i, j, nbNonZero;
    
    // Open file
    FILE* stream = fopen(path_adjancecy_matrix_file, "r");
    if(stream == NULL) {
        perror("Error opening file");
        exit(1);
    }

    // read header line;
    fgets(line, 1024, stream);
    sscanf(line, "%d %d %d", &nb_line, &nb_col, &nbNonZero);

    if(nb_line != nb_col){
        printf("Adjancecy matrix is not a square matrix:\nnb_line == %d != nb_col == %d\n", nb_line, nb_col);
        exit(1);
    }
    *matrix_dim = nb_line;

    // If sparce representation (compressed format only non zero element)
    if(sparce_rep){
        *sparce_adjacency_matrix = calloc(nbNonZero * 3, sizeof(int));
    }
    else{
        *adjacency_matrix = calloc(nb_line * nb_col, sizeof(unsigned char));
    }

    nbNonZero = 0;
    while (fgets(line, 1024, stream))
    {
        sscanf(line, "%d %d", &i, &j);

        if(sparce_rep){
            // insert element into the sparce matrix
            (*sparce_adjacency_matrix)[nbNonZero * 3] = i;
            (*sparce_adjacency_matrix)[nbNonZero * 3 + 1] = j;
            (*sparce_adjacency_matrix)[nbNonZero * 3 + 2] = 1;
            nbNonZero++;         
        }
        else{
            (*adjacency_matrix)[i * nb_line + j] = 1;
        }

    }

    fclose(stream);
    return nbNonZero;
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
double* create_transition_matrix(unsigned char* adjacency_matrix, int* sparce_adjacency_matrix,  int n, int nbNonZero, int sparce_rep, int parralel){
    
    double *out_links_vector = calloc(n,sizeof(double)),
            *vector_of_ones = calloc(n,sizeof(double)),
            *transition_matrix;
    int i;

    #pragma omp parallel for schedule(static) if(parralel)
    for(i = 0; i < n; i++){
        vector_of_ones[i] = 1;
    }

    // Get out_links_vectors
    if(sparce_rep){
        /*
            Time complexity : TO(nzero)
        */
        Sparce_Matrix_Vector_Product_int(sparce_adjacency_matrix, vector_of_ones, nbNonZero, n, out_links_vector, parralel);
        free(vector_of_ones);

        transition_matrix = calloc(nbNonZero * 3, sizeof(double));
        #pragma omp parallel for schedule(static) if(parralel)
        for(i = 0; i < nbNonZero; i++){
            // if link exists from j to i
            // then probabilty of going to i from j is 1 div number of out links of node j
            int j = sparce_adjacency_matrix[i*3];
            transition_matrix[i*3] = sparce_adjacency_matrix[i*3+1];
            transition_matrix[i*3+1] = j;
            transition_matrix[i*3+2] = 1 / out_links_vector[j];
        }
    }
    else{
        Matrix_Vector_Product_uchar(adjacency_matrix, vector_of_ones, n, n, out_links_vector, parralel);
        free(vector_of_ones);

        transition_matrix = calloc(n * n, sizeof(double));
        #pragma omp parallel for schedule(static) if(parralel)
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
        
    free(out_links_vector);

    return transition_matrix;
}