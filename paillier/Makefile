CC=gcc
CFLAGS=-Wall -Werror -ansi -pedantic -O3 -std=c11 -mrdrnd

all: paillier.o
	$(CC) $(CFLAGS) -o paillier paillier.o -lgmp

paillier.o: paillier.c
	$(CC) $(CFLAGS) -c paillier.c

clean:
	@(rm -f *.o)

remove:

	@(rm -f *.o paillier)
