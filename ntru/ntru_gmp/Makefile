CC=gcc
CFLAGS=-ansi -pedantic -Ofast -Wall -Werror -std=c11 -mrdrnd

# all: poly.o
	# $(CC) $(CFLAGS) -o poly poly.o -lgmp
all: ntru.o poly.o
	$(CC) $(CFLAGS) -o ntru ntru.o poly.o -lgmp

poly.o: poly.c poly.h
	$(CC) $(CFLAGS) -c poly.c

ntru.o: ntru.c ntru.h
	$(CC) $(CFLAGS) -c ntru.c

remove:
	@(rm -f *.o poly ntru)
