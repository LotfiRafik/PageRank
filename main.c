#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "hashtable.h"
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
    for (int i=0; i<size; i++)
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
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
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
    for (int i = 0; i < row; i++){
        for (int j = 0; j < col; j++){
            fprintf(stream,"%0.2f ", matrix[i*col+j]);
        }
        fprintf(stream,"\n");
    }

    fclose(stream);
}

// Return artist number, map artistId, indexID
int parse_infoArtist_csv(char* filepath, HashTable* hashtable){

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

        ht_insert(hashtable, artist_id, artist_cpt);

        artist_cpt++;
    }

    printf("Nombre d'artiste == %d\n", artist_cpt);

    // Close FIle
    fclose(stream);

    return 0;
}

// Return artist number, map artistId, indexID
int parse_collaborations_csv(char* filepath, HashTable* hashtable, double* adjacency_matrix, int nb_artists){

    FILE* stream = fopen(filepath, "r");
    char  delimiter = ';';
    int nb_csv_columns = 0;
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
        adjacency_matrix[source_idx * nb_artists +  target_idx] = 1;

    }

    // Close FIle
    fclose(stream);
}

double* create_transition_matrix(double* adjacency_matrix, double *out_links_vector, int nb_nodes){
    // Initialize transition matrix
    double *transition_matrix = calloc(nb_nodes * nb_nodes, sizeof(double));

    for (int i = 0; i < nb_nodes; i++){
        for (int j = 0; j < nb_nodes; j++){
            // if link exists from j to i
            // then probabilty of going to i from j is 1 div number of out links of node j
            if(adjacency_matrix[j*nb_nodes+i]){
                transition_matrix[i*nb_nodes+j] = 1 / out_links_vector[j];
            }

        }
    }

    return transition_matrix;
}

int main() {

    char* relative_path = "./Deezer-small-DS/";
    char datasets[][254] = {"Adele/", "Taylor Swift/", "David Guetta/", "Exo Td/"};
    int nb_artists_per_dataset[] = {27 , 220, 24851, 4};
    char files[][254] = {"InfoArtist.csv", "collaborations.csv", "level.csv"};

    // Dataset to use
    const int dataset_idx = 0;

    char file_path[254] = "";
    char artistInfo_file_path[254] = "";
    char collaborations_file_path[254] = "";
    char level_file_path[254] = "";

    strcat(file_path, relative_path);
    strcat(file_path, datasets[dataset_idx]);

    // copying str1 to str2
    strcpy(artistInfo_file_path, file_path);
    strcpy(collaborations_file_path, file_path);    
    strcpy(level_file_path, file_path);

    strcat(artistInfo_file_path, files[0]);
    strcat(collaborations_file_path, files[1]);
    strcat(level_file_path, files[2]);

    // Number of artists of the current used dataset
    int nb_artists = nb_artists_per_dataset[dataset_idx];

    // Create HashMap to map each artist ID to it virtuel index between [0, nb_artists]
    HashTable* hashtable = create_table(nb_artists);

    // Parse InfoArtist.csv file
    if(parse_infoArtist_csv(artistInfo_file_path, hashtable) < 0){
        perror("Error : can not parse InfoArtist file");
        exit(-1);
    }

    // Initialize adjacency matrix
    double *adjacency_matrix = calloc(nb_artists * nb_artists, sizeof(double));

    if(parse_collaborations_csv(collaborations_file_path, hashtable, adjacency_matrix, nb_artists) < 0){
        perror("Error : can not parse collaborations file");
        exit(-1);
    }

    char adjacency_matrix_file_path[254] = "";
    strcpy(adjacency_matrix_file_path, file_path);
    strcat(adjacency_matrix_file_path, "adjacency_matrix.csv");

    write_adjacency_matrix(nb_artists, nb_artists, adjacency_matrix, adjacency_matrix_file_path);

    double *out_links_vector = calloc(nb_artists,sizeof(double));
    double *vector_of_ones = calloc(nb_artists,sizeof(double));
    for (int i = 0; i < nb_artists; i++){
        vector_of_ones[i] = 1;
    }
    // Get out_links_vectors
    Matrix_Vector_Product(adjacency_matrix, vector_of_ones, nb_artists, nb_artists, out_links_vector);
    free(vector_of_ones);

    double* transition_matrix = create_transition_matrix(adjacency_matrix, out_links_vector, nb_artists);
    

    char transition_matrix_file_path[254] = "";
    strcpy(transition_matrix_file_path, file_path);
    strcat(transition_matrix_file_path, "transition_matrix.csv");

    write_matrix(nb_artists, nb_artists, transition_matrix, transition_matrix_file_path);

    // Free
    free(out_links_vector);
    free(adjacency_matrix);
    free(transition_matrix);
    // Free HashTable
    free_table(hashtable);
    return 0;
}