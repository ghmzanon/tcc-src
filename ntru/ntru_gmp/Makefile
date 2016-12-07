CC=gcc
CFLAGS=-ansi -pedantic -Ofast -Wall -Werror -std=c11 -mrdrnd

poly: poly_test.o poly.o
	$(CC) $(CFLAGS) -o poly_test poly_test.o poly.o -lgmp

ntru: ntru_test.o ntru.o poly.o
	$(CC) $(CFLAGS) -o ntru_test ntru_test.o ntru.o poly.o -lgmp

poly.o: poly.c poly.h
	$(CC) $(CFLAGS) -c poly.c

poly_test.o: poly_test.c
	$(CC) $(CFLAGS) -c poly_test.c

ntru.o: ntru.c ntru.h
	$(CC) $(CFLAGS) -c ntru.c

ntru_test.o: ntru_test.c
	$(CC) $(CFLAGS) -c ntru_test.c

remove:
	@(rm -f *.o poly_test ntru_test)
