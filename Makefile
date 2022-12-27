CC = gcc
EXEC = a.out

all: $(EXEC)

main.o: main.c
	$(CC) -c main.c

hashtable.o: hashtable.c
	$(CC) -c hashtable.c

blas.o: blas.c
	$(CC) -c blas.c

page_rank.o: page_rank.c
	$(CC) -c page_rank.c

$(EXEC): main.o hashtable.o blas.o page_rank.o
	$(CC) -o $(EXEC) main.o hashtable.o blas.o page_rank.o -fopenmp -lm

clean:
	rm *.o

cleanall:
	rm *.o
	rm $(EXEC)
