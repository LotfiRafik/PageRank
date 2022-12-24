# PageRank algorithm implemention in C

### Usage :  
Compile :  
gcc main.c hashtable.c blas.c page_rank.c -fopenmp -lm

Execute :  
./a.out dataset<0-3> parallel<0|1> sparce<0|1> outputFiles<0|1>

Datasets:  
| id |  name | artist's number |
| --- | ----------- | ---------- |
| 0 | Adele | 28 |
| 1 | Taylor Swift | 220 |
| 2 | David Guetta | 24851 |
| 3 | Exo Td | 4 |
