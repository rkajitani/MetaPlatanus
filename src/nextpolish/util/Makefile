CC = gcc
CFLAGS += -Wall -O3 -s
PROGS = seq_split seq_count

all: $(PROGS)

seq_split: kseq.h thpool.h seq_split.c thpool.c
	$(CC) $(CFLAGS) $(LDFLAGS) -pthread -o $@ seq_split.c thpool.c -lz 

seq_count: seq_count.c kseq.h
	$(CC) $(CFLAGS) $(LDFLAGS) -std=c99 -o $@ seq_count.c -lz

clean:
	rm -rf *.o $(PROGS)
