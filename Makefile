CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

all:lt-tag

lt-tag:lt-tag.c kthread.c kseq.h khash.h
		$(CC) $(CFLAGS) -pthread lt-tag.c kthread.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session* lt-tag
