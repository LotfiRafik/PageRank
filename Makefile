CC = gcc -Wno-unused-result -O3 
EXEC = a.out

all: $(EXEC)

deezer:  deezer.o utils/hashtable.o utils/blas.o utils/utils.o page_rank.o
	$(CC) -o $(EXEC) deezer.o utils/hashtable.o utils/blas.o utils/utils.o page_rank.o -fopenmp -lm

deezer.o: deezer.c
	$(CC) -c deezer.c -fopenmp

main.o: main.c
	$(CC) -c main.c -fopenmp

hashtable.o: utils/hashtable.c
	$(CC) -c hashtable.c -fopenmp

blas.o: utils/blas.c
	$(CC) -c blas.c -fopenmp

utils.o: utils/utils.c
	$(CC) -c utils/utils.c -fopenmp

page_rank.o: page_rank.c
	$(CC) -c page_rank.c -fopenmp

$(EXEC): main.o utils/hashtable.o utils/blas.o utils/utils.o page_rank.o
	$(CC) -o $(EXEC) main.o utils/hashtable.o utils/blas.o utils/utils.o page_rank.o -fopenmp -lm

clean:
	rm *.o

cleanall:
	rm *.o
	rm $(EXEC)
