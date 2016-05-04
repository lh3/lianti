CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
PROG=lt-trim

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lt-trim:lt-trim.o kthread.o kseq.h khash.h
		$(CC) $(CFLAGS) lt-trim.o kthread.o -o $@ -lz -lm -lpthread

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

lt-trim.o: kvec.h khash.h kseq.h
