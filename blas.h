double* initRandomVector(int n, int min, int max);

double* initRandomMatrix(double* A, int row,int col, int min, int max);

void displayVector(double* x, int n);

void displayMatrix(int row, int col, double* matrix);

double* transposeMatrix(int row, int col, double* A);

/*
    x= α.Ax+βy
    n*n : size matrix A
    nzero : size sparce matrix (i.e number of non-zero elements of matrix A)

    Time/Space complexity :
    Sequential case:
        using normal representation: TO(n^2) SO(n)
        using sparce representation: TO(nzero + n) SO(n)
    Parallel case:
        using normal representation: TO(n^2 / p) SO(n)
        using sparce representation: TO(nzero / p + n + n/p) SO(n*p + n)
*/
void blas21(double* A, double* x, double* y, double alpha, double beta, int row, int col, int nbNonZeroA, int parallel, int sparce_rep);

/*
    x= α.Ax+βy
    n : matrix size
    Time complexity : TO(n^2 / p)
    Space complexity : SO(n)
*/
void blas21_parallel(double* A, double* x, double* y, double alpha, double beta, int row, int col);
/*
    x= α.Ax+βy
    n : matrix size
    Time complexity : TO(n^2)
    Space complexity : SO(n)
*/
void blas21_sequential(double* A, double* x, double* y, double alpha, double beta, int row, int col);
/*
    x= α.Ax+βy
    p : nb processors
    n*n : size matrix
    nzero : size sparce matrix (i.e number of non-zero elements)
    Time complexity : TO(nzero / p + n + n/p)
    Space complexity : SO(n*p + n)
*/
void blas21_parallel_sparce(double* sparceA, double* x, double* y, double alpha, double beta, int n, int sizeSparceA);
/*
    x= α.Ax+βy
    n*n : size matrix
    nzero : size sparce matrix (i.e number of non-zero elements)
    Time complexity : TO(nzero + n)
    Space complexity : SO(n)
*/
void blas21_sequential_sparce(double* sparceA, double* x, double* y, double alpha, double beta, int n, int sizeSparceA);

void Vector_Scalar_Product(double* x, double alpha, int n);
void Vector_Scalar_Product_parallel(double* x, double alpha, int n);

void Vector_Vector_Addition(double* x, double* y, int n);

/* 
    Produit scalaire 2 vecteurs
    p : nb processors
    Time/Space complexity :
        Sequential case: TO(n) SO(1)
        Parallel case: TO(n / p) SO(1)
*/
double DotProduct(double* x, double* y, int n, int parallel);

void Matrix_Vector_Product(double* A, double* v, int row, int col, double* Av, int parallel);

void Matrix_Vector_Product_sequential(double* A, double* v, int row, int col, double* Av);
void Matrix_Vector_Product_parralel(double* A, double* v, int row, int col, double* Av);

double* Matrix_Matrix_Product(double* A, double* B, int rowA, int colA, int colB, double* AB);
double* Matrix_Matrix_Subsctraction(double* A, double* B, int row, int col, double* AB);

/*
    Norme2 d'un vecteur
    p : nb processors
    Time/Space complexity :
        Sequential case: TO(n) SO(1)
        Parallel case: TO(n / p) SO(1)
*/
double Norme(double* x, int n, int parallel);



// La norme Frobenius d’une matrice
double NormeFrobenius(int nb_ligne, int nb_colonne, double* matrice);




/*   
    n*n : size square matrix
    nzero : size sparce matrix (i.e number of non-zero elements in the original matrix)
    Sequential case:
        Time complexity :  TO(nzero)
    Parallel case:
        Time complexity : TO(nzero / p + n*p)
        Space complexity : SO(n * p)
*/
void Sparce_Matrix_Vector_Product(double* A, double* v, int sizeA, int sizeV, double* Av, int parallel);

int count_zero_matrix(double *A, int row, int col);

/* 
    Convert normal matrix representation to sparce matrix representation
    Time/Space complexity:
        sequential case: TO(n^2), SO(nzero * 3)
        parallel case:
*/
double* matrix_to_sparce(double *A, int row, int col, int* nbNonZero);

// Convert sparce matrix representation to normal matrix representation
double* sparce_to_matrix(double *sparceA, double *A, int row, int col, int nbNonZero);

void display_sparce_matrix(double *sparceA, int nbNonZero);
