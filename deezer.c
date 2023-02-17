#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include <omp.h>
#include "utils/hashtable.h"
#include "utils/blas.h"
#include "utils/utils.h"
#include "page_rank.h"

// *************************** Functions Declarations *****************************
int parse_infoArtist_csv(char* filepath, HashTable** hashtable, double** fans_percent_artists, int** artists_ids, char*** artists_names);
int number_collaborations_csv(char* filepath, HashTable* hashtable);
int parse_collaborations_csv(char* filepath, HashTable* hashtable, unsigned char** adjacency_matrix, int** sparce_adjacency_matrix, int nb_artists, int sparce_rep);
// ***************************************************************************************


// ***************************** MAIN *****************************
int main(int argc, char **argv) {
    double temps, debut,fin,
            *fans_percent_artists,
            *opposite_rep_adjacency_matrix;
    unsigned char *adjacency_matrix;
    int *sparce_adjacency_matrix,
        i,
        *artists_ids;

    HashTable* hashtable;

    char *relative_path = "./Deezer-DS/",
        datasets[][254] = {"Adele/", "Taylor Swift/", "David Guetta/", "Ed Sheeran/"},
        files[][254] = {"InfoArtist.csv", "collaborations.csv", "level.csv"},
        file_path[254] = "",
        artistInfo_file_path[254] = "",
        collaborations_file_path[254] = "",
        level_file_path[254] = "",
        **artists_names;

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
    strcpy(artistInfo_file_path, file_path);
    strcpy(collaborations_file_path, file_path);    
    strcpy(level_file_path, file_path);
    strcat(artistInfo_file_path, files[0]);
    strcat(collaborations_file_path, files[1]);
    strcat(level_file_path, files[2]);


    // ************************ BEGIN : Parse InfoArtist.csv file *****************************
    debut = omp_get_wtime();

    int nb_artists = parse_infoArtist_csv(artistInfo_file_path, &hashtable, &fans_percent_artists, &artists_ids, &artists_names);
    if(nb_artists < 0){
        perror("Error : can not parse InfoArtist file");
        exit(-1);
    }

    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to parse infoArtist.csv ==>  : %lf s\n",temps);
    // ************************ END : Parse InfoArtist.csv file *****************************


    // ************************ BEGIN : Parse Collaborations.csv file *****************************
	debut = omp_get_wtime();

    // Initialize adjacency matrix 
    int nb_collaborations = parse_collaborations_csv(collaborations_file_path, hashtable, &adjacency_matrix, &sparce_adjacency_matrix, nb_artists, SPARCE_REPRESENTATION);
    free_table(hashtable);
    if(nb_collaborations < 0){
        perror("Error : can not parse collaborations file");
        exit(-1);
    }

    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create adjacency_matrix ==>  : %lf s\n",temps);
    // ************************ END : Parse Collaborations.csv file *****************************


    // ************************ BEGIN : CREATE TRANSITION MATRIX **************************************
    debut = omp_get_wtime();
    
    double* transition_matrix = create_transition_matrix(adjacency_matrix, sparce_adjacency_matrix, nb_artists, nb_collaborations, SPARCE_REPRESENTATION, MODE_EXEC);
    
    fin = omp_get_wtime();
	temps = fin - debut;
	printf("Time to create transition_matrix ==>  : %lf s\n",temps);
    // ************************ END : CREATE TRANSITION MATRIX **************************************


    // ************************ BEGIN : Apply PageRank Algorithm **************************************
	debut = omp_get_wtime();
    // {
    //     #pragma omp parallel for schedule(static) private(i)
    //     for (i = 0; i < nb_artists; i++)
    //     {
    //         // initializer le vecteur de depart
    //         fans_percent_artists[i] = 1./nb_artists;
    //     }    
    // }
    double *pg_vector = page_rank(transition_matrix, fans_percent_artists, nb_collaborations, nb_artists, 0.85, 0.00001, MODE_EXEC, SPARCE_REPRESENTATION);
	fin = omp_get_wtime();
	temps = fin - debut;

	printf("Time PageRank Algorithm  ==>  : %lf s\n",temps);
    // ************************ END : Apply PageRank Algorithm **************************************



    // ************************ BEGIN : Write Results to DISK ************************

    // Ouput pagerank score of each artist
    char artists_scores_file_path[254] = "";
    strcpy(artists_scores_file_path, file_path);
    strcat(artists_scores_file_path, "artists_scores.csv");
    FILE* stream = fopen(artists_scores_file_path, "w");
    if(stream == NULL) {
        perror("Error opening artists_scores file");
        exit(1);
    }
    fprintf(stream,"artist_id;artist_name;pourtentage_fans;pagerank_score\n");
    for (size_t i = 0; i < nb_artists; i++){      
        fprintf(stream,"%d;%s;%0.2f;%0.15f\n", artists_ids[i], artists_names[i], fans_percent_artists[i]*100 ,pg_vector[i]);
    }
    fclose(stream);

    if(outputfiles){
        // WRITE ADJACENCY MATRIX TO DISK
        char adjacency_matrix_file_path[254] = "";
        strcpy(adjacency_matrix_file_path, file_path);
        strcat(adjacency_matrix_file_path, "adjacency_matrix.csv");
        write_adjacency_matrix(nb_artists, nb_artists, nb_collaborations, adjacency_matrix, sparce_adjacency_matrix, adjacency_matrix_file_path, SPARCE_REPRESENTATION);      

        // WRITE TRANSITION MATRIX TO DISK
        char transition_matrix_file_path[254] = "";
        strcpy(transition_matrix_file_path, file_path);
        strcat(transition_matrix_file_path, "transition_matrix.csv");
        write_matrix(nb_artists, nb_artists, nb_collaborations, transition_matrix, transition_matrix_file_path, SPARCE_REPRESENTATION);
    }


    // ALWAYS WRITE PAGERANK VECTOR TO DISK
    char pagerank_vector_file_path[254] = "";
    strcpy(pagerank_vector_file_path, file_path);
    strcat(pagerank_vector_file_path, "pagerank_vector.csv");
    write_vector(nb_artists, pg_vector, pagerank_vector_file_path);


    // ************************ END : Write Results to DISK ************************


    // Free allocated space
    if(SPARCE_REPRESENTATION)
        free(sparce_adjacency_matrix);
    else
        free(adjacency_matrix);

    free(transition_matrix);

    return 0;
}

/*  
    map each artistId to an integer between (0, nb_artists)
    returns integer : if >= 0 : nb artists (i.e nb lines in the csv file) OR < 0 : error 
*/
int parse_infoArtist_csv(char* filepath, HashTable** hashtable, double** fans_percent_artists, int** artists_ids, char*** artists_names){

    FILE* stream = fopen(filepath, "r");
    char delimiter = ';';
    int nb_csv_columns = 0,
        total_nb_fans = 0,
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
        int nb_fans = atoi(tokens[3]);
        if(!artist_id){
            perror("Can not convert artist id to integer");
            return(-1);
        }
        
        free_splitted_tokens(tokens, nb_csv_columns);

        total_nb_fans += nb_fans;
        artist_cpt++;
    }

    printf("Nombre d'artiste == %d\n", artist_cpt);

    // Create HashMap to map each artist ID to it virtuel index between [0, nb_artists]
    *hashtable = create_table(artist_cpt);

    // Create vector of % of fans of each artist
    *fans_percent_artists = malloc(artist_cpt * sizeof(double));
    
    // Create vector of artists ids
    *artists_ids = malloc(artist_cpt * sizeof(int));
    *artists_names = malloc(artist_cpt * sizeof(char*));

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
        int nb_fans = atoi(tokens[3]);
        if(!artist_id){
            perror("Can not convert artist id to integer");
            return(-1);
        }
        
        // map artist_id to index
        
        ht_insert(*hashtable, (long)artist_id, artist_cpt);
        // map index to artist_id
        (*artists_ids)[artist_cpt] = artist_id;

        // Vector of artist names
        (*artists_names)[artist_cpt] = strdup(tokens[2]);
        
        (*fans_percent_artists)[artist_cpt] = (double) nb_fans / (double) total_nb_fans;
        artist_cpt++;

        free_splitted_tokens(tokens, nb_csv_columns);
    }


    // Close FIle
    fclose(stream);

    return artist_cpt;
}

/*
    return number of artists's collaborations
*/
int number_collaborations_csv(char* filepath, HashTable* hashtable){

    FILE* stream = fopen(filepath, "r");
    char  delimiter = ';';
    int nb_csv_columns = 0;
    int nb_collaborations = 0;
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
        int source_idx = ht_search(hashtable, (long)source_id);
        int target_idx = ht_search(hashtable, (long)target_id);

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


        nb_collaborations++;     
    }

    // Close FIle
    fclose(stream);

    return nb_collaborations;
}

/*
    create adjacency matrix from collaborations csv file
    return nb of collaborations (nb non zero elements of the adjacency matrix)
*/
int parse_collaborations_csv(char* filepath, HashTable* hashtable, unsigned char** adjacency_matrix, int** sparce_adjacency_matrix, int nb_artists, int sparce_rep){

    char  delimiter = ';',
                line[1024];
    int nb_csv_columns = 0,
        nb_not_zero_adjacency_matrix = 0,
        line_cpt = 0,
        nb_collaborations = -1;
    
    FILE* stream = fopen(filepath, "r");
    if(stream == NULL) {
        perror("Error opening file");
        return(-1);
    }

    // @FIX could contain duplicate collaborations
    nb_collaborations = number_collaborations_csv(filepath, hashtable);

    // Create HashMap to save collaborations
    HashTable *collab_hashtable = create_table(nb_collaborations);

    // If sparce representation (compressed format i.e only non zero element)
    if(sparce_rep)
        *sparce_adjacency_matrix = calloc(nb_collaborations * 3, sizeof(int));
    else
        *adjacency_matrix = calloc(nb_artists * nb_artists, sizeof(unsigned char));
    
    // reset nb_collaborations
    nb_collaborations = 0;

    // read header line;
    fgets(line, 1024, stream);
    // parse line by line
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
        int source_idx = ht_search(hashtable, (long)source_id);
        int target_idx = ht_search(hashtable, (long)target_id);

        if(source_idx < 0){
            // printf("Key %d not found in hashtable", source_id);
            // return(-1);
            continue;
        }
        if(target_idx < 0){
            // printf("Key %d not found in hashtable", target_id);
            // return(-1);
            // Skip unfound artist
            continue;
        }

        // Check if we already inserted this collaboration
        // compute unique integer from two integers http://en.wikipedia.org/wiki/pairing_function

        long s = (long)source_idx;
        long t = (long)target_idx;
        // unsigned long unique_id = (source_idx + target_idx)*(source_idx + target_idx + 1)/2 + target_idx;
        long unique_id = (s + t)*(s + t + 1)/2 + t;
        if(unique_id < 0){
            continue;
        }
        if(ht_search(collab_hashtable, unique_id) == -1){
            // insert collaboration into hashmap
            ht_insert(collab_hashtable, unique_id, 1); 
            if(sparce_rep){
                // insert element into the sparce matrix
                (*sparce_adjacency_matrix)[nb_collaborations * 3] = source_idx;
                (*sparce_adjacency_matrix)[nb_collaborations * 3 + 1] = target_idx;
                (*sparce_adjacency_matrix)[nb_collaborations * 3 + 2] = 1;   
                
            }else
                (*adjacency_matrix)[source_idx * nb_artists +  target_idx] = 1;
            
            nb_collaborations++;
        }

        // if((*adjacency_matrix)[source_idx * nb_artists +  target_idx] != 1){
        //     (*adjacency_matrix)[source_idx * nb_artists +  target_idx] = 1;
            
        //     if(sparce_rep){
        //         // insert element into the sparce matrix
        //         (*sparce_adjacency_matrix)[nb_not_zero_adjacency_matrix * 3] = source_idx;
        //         (*sparce_adjacency_matrix)[nb_not_zero_adjacency_matrix * 3 + 1] = target_idx;
        //         (*sparce_adjacency_matrix)[nb_not_zero_adjacency_matrix * 3 + 2] = 1;            
        //     }

        //     // Increment cpt non zero elements 
        //     // usefull to create sparce matrix representation later
        //     nb_not_zero_adjacency_matrix++;
        // }
    }

    // Close FIle
    fclose(stream);

    free_table(collab_hashtable);

    return nb_collaborations;
}