CC=gcc

make: 
	$(CC) QuicKmer.c -O3 -pthread -std=c99 -lm -o quicKmer
