char **split(char *string, char seperator, int* nb_tokens);
void free_splitted_tokens(char** tokens, int size);
void write_adjacency_matrix(int row, int col, int nbNonZero, unsigned char* adjacency_matrix, int* sparce_adjacency_matrix, char* filepath, int sparce_rep);
void write_matrix(int row, int col, int nbNonZero, double* matrix, char* filepath, int sparce_rep);
void write_vector(int dim, double* vector, char* filepath);
double* create_transition_matrix(unsigned char* adjacency_matrix, int* sparce_adjacency_matrix, int n, int nbNonZero, int sparce_rep, int parralel);
/*
    load adjacency matrix from file
    return nb non zero elements of the adjacency matrix
*/
int read_adj_matrix_from_file(char* path_adjancecy_matrix_file, unsigned char** adjacency_matrix, int** sparce_adjacency_matrix, int *matrix_dim, int sparce_rep);
