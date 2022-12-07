double* initRandomVector(int n, int min, int max);

double* initRandomMatrix(double* A, int row,int col, int min, int max);

void displayVector(double* x, int n);

void displayMatrix(int row, int col, double* matrix);

double* transposeMatrix(int row, int col, double* A);

void blas21(double* A, double* x, double* y, double alpha, double beta, int row, int col, int parallel);
void blas21_parallel(double* A, double* x, double* y, double alpha, double beta, int row, int col);
void blas21_sequential(double* A, double* x, double* y, double alpha, double beta, int row, int col);

void Vector_Scalar_Product(double* x, double alpha, int n);
void Vector_Scalar_Product_parallel(double* x, double alpha, int n);

void Vector_Vector_Addition(double* x, double* y, int n);

// Produit scalaire 2 vecteurs
double DotProduct(double* x, double* y, int n);

double DotProduct_parallel(double* x, double* y, int n);

void Matrix_Vector_Product(double* A, double* v, int row, int col, double* Av, int parallel);

void Matrix_Vector_Product_sequential(double* A, double* v, int row, int col, double* Av);
void Matrix_Vector_Product_parralel(double* A, double* v, int row, int col, double* Av);

double* Matrix_Matrix_Product(double* A, double* B, int rowA, int colA, int colB, double* AB);
double* Matrix_Matrix_Subsctraction(double* A, double* B, int row, int col, double* AB);
// Norme2 d'un vector
double Norme(double* x, int n);



// La norme Frobenius d’une matrice
double NormeFrobenius(int nb_ligne, int nb_colonne, double* matrice);