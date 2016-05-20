CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
PROG=lianti
OBJS=kthread.o bgzf.o hts.o bedidx.o sam.o \
	 trim.o ldup.o group.o count.o cnv.o \
	 main.o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lianti:$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $@ -lz -lm -lpthread

clean:
		rm -fr gmon.out *.o ext/*.o a.out *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bedidx.o: ksort.h kseq.h khash.h
bgzf.o: bgzf.h
cnv.o: kvec.h kseq.h ksort.h
count.o: kvec.h kseq.h kdq.h
group.o: sam.h bgzf.h hts.h kdq.h kvec.h ksort.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
ldup.o: sam.h bgzf.h hts.h kdq.h khash.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h
trim.o: kvec.h khash.h kseq.h
