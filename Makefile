CC=gcc
CFLAGS=-I.
DEPS=hashmap.h
OBJ=hashmap.o histo-hash.o

all: histo-hash histo-vector

%.o: %.c $(DEPS)
	$(CC) -Wall -c -o $@ $< $(CFLAGS)

histo-hash: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

histo-vector: histo-vector.c
	$(CC) -o $@ $^ $(CFLAGS) -lm

clean:
	rm -f histo-vector histo-hash  $(OBJ) *~ 

