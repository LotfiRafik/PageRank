#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "hashtable.h"
#include "blas.h"

char** parse_csv_line(char* line, int nb_columns, char* delimiter){

    // Init list of tokens/strings
    char** tokens = calloc(nb_columns, sizeof(char*));
    // Init nb token 
    int counter = 0;

    char* strToken = strtok(line, delimiter);

    while ( strToken != NULL && counter < nb_columns) {
        tokens[counter] = strdup(strToken);
        // Get next token.
        strToken = strtok(NULL, delimiter);
        counter++;
    }

    return tokens;
}

void free_csv_line(char** tokens, int size){
    for (int i=0; i<size; i++){
        if (tokens[i] != NULL){
            free(tokens[i]);
        }
    }
    free(tokens);
}


// Return artist number, map artistId, indexID
int parse_infoArtist_csv(char* filepath, HashTable* hashtable){

    FILE* stream = fopen(filepath, "r");
    char* delimiter = ";";
    int nb_csv_columns = 5;
    int artist_cpt = 0;
    if(stream == NULL) {
        perror("Error opening file");
        return(-1);
    }

    char line[1024];
    // read header line;
    fgets(line, 1024, stream);
    while (fgets(line, 1024, stream))
    {   
        char**  tokens = parse_csv_line(line, nb_csv_columns, delimiter);
        // for (int i=0; i<nb_csv_columns; i++)
        //     printf("%s \n",tokens[i]);
        
        int artist_id = atoi(tokens[1]);
        if(!artist_id){
            perror("Can not convert artist id to integer");
            return(-1);
        }
        // printf("%s \n",tokens[1]);
        ht_insert(hashtable, artist_id, artist_cpt);

        artist_cpt++;
        free_csv_line(tokens, nb_csv_columns);
    }

    printf("Nombre d'artiste == %d\n", artist_cpt);

    // Close FIle
    fclose(stream);

    return 0;
}

// Return artist number, map artistId, indexID
int parse_collaborations_csv(char* filepath, HashTable* hashtable, double* adjacency_matrix, int nb_artists){

    FILE* stream = fopen(filepath, "r");
    char* delimiter = ";";
    int nb_csv_columns = 8;
    if(stream == NULL) {
        perror("Error opening file");
        return(-1);
    }

    char line[1024];
    // read header line;
    fgets(line, 1024, stream);
    while (fgets(line, 1024, stream))
    {   
        char**  tokens = parse_csv_line(line, nb_csv_columns, delimiter);
        int source_id = atoi(tokens[4]);
        int target_id = atoi(tokens[6]);
        if(!source_id || !target_id){
            perror("Can not convert source or target ID to integer");
            return(-1);
        }

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

        free_csv_line(tokens, nb_csv_columns);
    }

    // Close FIle
    fclose(stream);
}


int main() {

    char* relative_path = "./Deezer-small-DS/";
    char datasets[][254] = {"Adele/", "Taylor Swift/", "David Guetta/"};
    int nb_artists_per_dataset[] = {27 , 220, 24851};
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


    displayMatrix(nb_artists, nb_artists, adjacency_matrix);

    // Free adje
    free(adjacency_matrix);
    // Free HashTable
    free_table(hashtable);
    return 0;
}