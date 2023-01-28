#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "hashtable.h"
#include "blas.h"
#include "page_rank.h"
#include <omp.h>

// ******* Declarations *******

char **split(char *string, char seperator, int* nb_tokens);
void free_splitted_tokens(char** tokens, int size);
void write_adjacency_matrix(int row, int col, double* matrix, char* filepath);
void write_matrix(int row, int col, double* matrix, char* filepath);
int parse_infoArtist_csv(char* filepath, HashTable** hashtable);
int parse_collaborations_csv(char* filepath, HashTable* hashtable, double* adjacency_matrix, int nb_artists);
double* create_transition_matrix(double* adjacency_matrix, double *out_links_vector, int n, int nbNonZero, int sparce_rep, int parralel);

// *****************************


// ******* MAIN *******
int main(int argc, char **argv) {
    double temps, debut,fin,
            *opposite_rep_adjacency_matrix;
    HashTable* hashtable;

    char* relative_path = "./Deezer-small-DS/";
    char datasets[][254] = {"Adele/", "Taylor Swift/", "David Guetta/", "Ed Sheeran/", "Exo Td/"};
    char files[][254] = {"InfoArtist.csv", "collaborations.csv", "level.csv"};

    char file_path[254] = "";
    char artistInfo_file_path[254] = "";
    char collaborations_file_path[254] = "";
    char level_file_path[254] = "";



    // Dataset to use
    int dataset_idx = 0;
    // 1 : Parallel , 0 : Sequential
    int MODE_EXEC = 1;
    // 1 : Use sparce matrix representation, 0 : Use normal matrix representation
    int SPARCE_REPRESENTATION = 1;

    if(argc < 5){
        fprintf(stderr, "Usage: %s dataset<0-3> parallel<0|1> sparce<0|1> outputFiles<0|1>\n", argv[0]);
        exit(1);
    }
    dataset_idx = atoi(argv[1]);
    MODE_EXEC = atoi(argv[2]);
    SPARCE_REPRESENTATION = atoi(argv[3]);
    int outputfiles = atoi(argv[4]);


    // Prepare file paths
    strcat(file_path, relative_path);
    strcat(file_path, datasets[dataset_idx]);
    // copying str1 to str2
    strcpy(artistInfo_file_path, file_path);
    strcpy(collaborations_file_path, file_path);    
    strcpy(level_file_path, file_path);
    strcat(artistInfo_file_path, files[0]);
    strcat(collaborations_file_path, files[1]);
    strcat(level_file_path, files[2]);


	debut = omp_get_wtime();
    // Parse InfoArtist.csv file
    int nb_artists = parse_infoArtist_csv(artistInfo_file_path, &hashtable);
    if(nb_artists < 0){
        perror("Error : can not parse InfoArtist file");
        exit(-1);
    }
    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create HashTable ==>  : %lf s\n",temps);


	debut = omp_get_wtime();
    // Initialize adjacency matrix 
    // SO(n^2)
    double *adjacency_matrix = calloc(nb_artists * nb_artists, sizeof(double));
    int nb_collaborations = parse_collaborations_csv(collaborations_file_path, hashtable, adjacency_matrix, nb_artists);
    if(nb_collaborations < 0){
        perror("Error : can not parse collaborations file");
        exit(-1);
    }

    // Create sparce matrix representation of the adjacency matrix
    if(SPARCE_REPRESENTATION){
        opposite_rep_adjacency_matrix = adjacency_matrix;
        // TO(n^2), SO(nzero * 3)
        adjacency_matrix = matrix_to_sparce(adjacency_matrix, nb_artists, nb_artists, &nb_collaborations);
    }


    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create adjacency_matrix ==>  : %lf s\n",temps);


    // CREATE TRANSITION MATRIX
    debut = omp_get_wtime();
    // SO(n)
    double *out_links_vector = calloc(nb_artists,sizeof(double));
    // SO(n)
    double *vector_of_ones = calloc(nb_artists,sizeof(double));

    int i;

    // TO(n / p)
    #pragma omp parallel for schedule(static) if(MODE_EXEC)
    for(i = 0; i < nb_artists; i++){
        vector_of_ones[i] = 1;
    }

    // Get out_links_vectors
    if(SPARCE_REPRESENTATION)
        /*
            Time complexity : TO(nzero)
        */
        Sparce_Matrix_Vector_Product(adjacency_matrix, vector_of_ones, nb_collaborations, nb_artists, out_links_vector, MODE_EXEC);
    else
        /*
            Sequential case: TO(n^2), SO(1)
            Parallel case: TO(n^2 / p), SO(1)
        */
        Matrix_Vector_Product(adjacency_matrix, vector_of_ones, nb_artists, nb_artists, out_links_vector, MODE_EXEC);
        
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
    double* transition_matrix = create_transition_matrix(adjacency_matrix, out_links_vector, nb_artists, nb_collaborations, SPARCE_REPRESENTATION, MODE_EXEC);
    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create transition_matrix ==>  : %lf s\n",temps);


    // Apply PageRank Algorithm
    /*
        Time/Space complexity :
        Sequential case:
            Normal representation: TO(nb_iterations * n^2) SO(4n)
            Sparce representation: TO(nb_iterations * nzero) SO(4n)
        Parallel case:
            Normal representation: TO(nb_iterations * (n^2 / p)) SO(4n)
            Sparce representation: TO(nb_iterations * (nzero + n/p)) SO(4n)
    */
	debut = omp_get_wtime();
    double *pg_vector = page_rank(transition_matrix, nb_collaborations, nb_artists, 0.85, 0.00001, MODE_EXEC, SPARCE_REPRESENTATION);
	fin = omp_get_wtime();
	temps = fin - debut;

	printf("Time PageRank Algorithm  ==>  : %lf s\n",temps);


    // WRITE MATRIX TO DISK
    if(outputfiles){
        char adjacency_matrix_file_path[254] = "";
        strcpy(adjacency_matrix_file_path, file_path);
        strcat(adjacency_matrix_file_path, "adjacency_matrix.csv");
        if(SPARCE_REPRESENTATION)
            write_adjacency_matrix(nb_artists, nb_artists, opposite_rep_adjacency_matrix, adjacency_matrix_file_path);
        else
            write_adjacency_matrix(nb_artists, nb_artists, adjacency_matrix, adjacency_matrix_file_path);      

        // WRITE TRANSITION MATRIX TO DISK
        char transition_matrix_file_path[254] = "";
        strcpy(transition_matrix_file_path, file_path);
        strcat(transition_matrix_file_path, "transition_matrix.csv");

        if(SPARCE_REPRESENTATION)
            write_matrix(nb_artists, nb_artists, sparce_to_matrix(transition_matrix, NULL, nb_artists, nb_artists, nb_collaborations), transition_matrix_file_path);
        else
            write_matrix(nb_artists, nb_artists, transition_matrix, transition_matrix_file_path);
    }


    // WRITE PAGERANK VECTOR TO DISK
    char pagerank_vector_file_path[254] = "";
    strcpy(pagerank_vector_file_path, file_path);
    strcat(pagerank_vector_file_path, "pagerank_vector.csv");
    write_matrix(1, nb_artists, pg_vector, pagerank_vector_file_path);


    // Free allocated space
    free(out_links_vector);
    free(adjacency_matrix);
    free(transition_matrix);
    free_table(hashtable);

    return 0;
}



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

// Write adjacency matrix to a file
void write_adjacency_matrix(int row, int col, double* matrix, char* filepath){
    FILE* stream = fopen(filepath, "w");
    if(stream == NULL) {
        perror("Error opening file");
        exit(1);
    }
    int i;
for(i = 0; i < row; i++){
        int j;
for(j = 0; j < col; j++){
            fprintf(stream,"%d ", (int) matrix[i*col+j]);
        }
        fprintf(stream,"\n");
    }

    fclose(stream);
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
    map each artistId to an integer between (0, nb_artists)
    returns integer : if >= 0 : nb artists (i.e nb lines in the csv file) OR < 0 : error 
*/
int parse_infoArtist_csv(char* filepath, HashTable** hashtable){

    FILE* stream = fopen(filepath, "r");
    char delimiter = ';';
    int nb_csv_columns = 0,
        artist_cpt = 0;
    char** tokens;
    if(stream == NULL) {
        perror("Error opening file");
        return(-1);
    }

    char line[1024];

    // Count nb lines/nb artists
    // read header line;
    fgets(line, 1024, stream);
    while (fgets(line, 1024, stream))
    {   
        tokens = split(line, delimiter, &nb_csv_columns);
        int artist_id = atoi(tokens[1]);
        if(!artist_id){
            perror("Can not convert artist id to integer");
            return(-1);
        }
        
        free_splitted_tokens(tokens, nb_csv_columns);

        artist_cpt++;
    }

    printf("Nombre d'artiste == %d\n", artist_cpt);

    // Create HashMap to map each artist ID to it virtuel index between [0, nb_artists]
    *hashtable = create_table(artist_cpt);

    // Rewind stream to the begining
    rewind(stream);
    artist_cpt = 0;

    // Parse file

    // read header line;
    fgets(line, 1024, stream);
    while (fgets(line, 1024, stream))
    {   
        tokens = split(line, delimiter, &nb_csv_columns);
        int artist_id = atoi(tokens[1]);
        if(!artist_id){
            perror("Can not convert artist id to integer");
            return(-1);
        }
        
        free_splitted_tokens(tokens, nb_csv_columns);

        ht_insert(*hashtable, artist_id, artist_cpt);

        artist_cpt++;
    }


    // Close FIle
    fclose(stream);

    return artist_cpt;
}

/*
    create adjacency matrix from collaborations csv file
    return nb of collaborations (nb non zero elements of the adjacency matrix)
*/
int parse_collaborations_csv(char* filepath, HashTable* hashtable, double* adjacency_matrix, int nb_artists){

    FILE* stream = fopen(filepath, "r");
    char  delimiter = ';';
    int nb_csv_columns = 0;
    int nb_not_zero_adjacency_matrix = 0;
    if(stream == NULL) {
        perror("Error opening file");
        return(-1);
    }

    char line[1024];
    int line_cpt = 0;
    // read header line;
    fgets(line, 1024, stream);
    while (fgets(line, 1024, stream))
    {   
        line_cpt++;
        char** tokens = split(line, delimiter, &nb_csv_columns);
        int source_id = atoi(tokens[4]);
        int target_id = atoi(tokens[6]);
        if(!source_id || !target_id){
            printf("Can not convert sourceID : %s or targetID : %s to integer in line %d\n", tokens[4], tokens[6], line_cpt);
            return(-1);
        }
        
        free_splitted_tokens(tokens, nb_csv_columns);

        // Get index from artistId
        int source_idx = ht_search(hashtable, source_id);
        int target_idx = ht_search(hashtable, target_id);

        if(source_idx == -1){
            // printf("Key %d not found in hashtable", source_id);
            // return(-1);
            continue;
        }
        if(target_idx ==  -1){
            // printf("Key %d not found in hashtable", target_id);
            // return(-1);
            // Skip unfound artist
            continue;
        }

        // Source points to Target
        if(adjacency_matrix[source_idx * nb_artists +  target_idx] != 1){
            adjacency_matrix[source_idx * nb_artists +  target_idx] = 1;
            // Increment cpt non zero elements 
            // usefull to create sparce matrix representation later
            nb_not_zero_adjacency_matrix++;
        }

    }

    // Close FIle
    fclose(stream);

    return nb_not_zero_adjacency_matrix;
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
    int i;

    if(sparce_rep){
        // SO(nbNonZero * 3)
        transition_matrix = calloc(nbNonZero * 3, sizeof(double));
        // TO(nzero / p)
        #pragma omp parallel for schedule(static) if(parralel)
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
