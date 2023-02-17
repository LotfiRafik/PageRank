# PageRank algorithm implemention in C

### Usage :  
Compile :  

#### To apply pageRank algo on an adjacency matrix from a file :

Compile:  
make   
Execute :   
./a.out path_adjancecy_matrix_file parallel<0|1> sparce<0|1> outputFiles<0|1>

example:  
./a.out data/CQ_MATRICES/100_diag_0.150000_perturbed_7500.mtx 1 1 0

#### To apply pageRank algo on deezer dataset :
Compile:  
make deezer  
Execute :    
./a.out dataset<0-3> parallel<0|1> sparce<0|1> outputFiles<0|1>

example:  
./a.out 2 1 1 0

Deezer Datasets:  
| id |  name | artist's number |
| --- | ----------- | ---------- |
| 0 | Adele | 28 |
| 1 | Taylor Swift | 220 |
| 2 | David Guetta | 24851 |
| 3 | Ed Sheeran | 502317 |
