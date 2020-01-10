CC=gcc

make: 
	$(CC) QuicKmer.c -O3 -g -pthread -std=c99 -lm -o quicKmer2
