CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
PROG=lt-trim lt-ldup lt-group

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lt-trim:lt-trim.o kthread.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

lt-ldup:lt-ldup.o bgzf.o hts.o sam.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

lt-group:lt-group.o bgzf.o hts.o sam.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

clean:
		rm -fr gmon.out *.o ext/*.o a.out *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bgzf.o: bgzf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
lt-group.o: sam.h bgzf.h hts.h kdq.h
lt-ldup.o: sam.h bgzf.h hts.h kdq.h khash.h
lt-trim.o: kvec.h khash.h kseq.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h
