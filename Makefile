CC = gcc
EXEC = a.out

all: $(EXEC)

perf:  perf.o hashtable.o blas.o page_rank.o
	$(CC) -o perf perf.o hashtable.o blas.o page_rank.o -fopenmp -lm

perf.o: perf.c
	$(CC) -c perf.c -fopenmp

main.o: main.c
	$(CC) -c main.c -fopenmp

hashtable.o: hashtable.c
	$(CC) -c hashtable.c -fopenmp

blas.o: blas.c
	$(CC) -c blas.c -fopenmp

page_rank.o: page_rank.c
	$(CC) -c page_rank.c -fopenmp

$(EXEC): main.o hashtable.o blas.o page_rank.o
	$(CC) -o $(EXEC) main.o hashtable.o blas.o page_rank.o -fopenmp -lm

clean:
	rm *.o

cleanall:
	rm *.o
	rm $(EXEC)
