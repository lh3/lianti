#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

typedef struct {
	int ctg;
	uint32_t lo:1, st:31;
	uint32_t ro:1, en:31;
	int n_seg, n_frag;
} lt_frag_t;

#include "kdq.h"
KDQ_INIT(lt_frag_t)

typedef struct {
	kstream_t *ks;
	gzFile fp;
	kstring_t s;
	int n_ctg, m_ctg;
	char **ctg;
} lt_reader_t;

lt_reader_t *lt_cnt_open(const char *fn)
{
	lt_reader_t *r;
	gzFile fp;
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	r = (lt_reader_t*)calloc(1, sizeof(lt_reader_t));
	r->fp = fp;
	r->ks = ks_init(fp);
	return r;
}

void lt_cnt_close(lt_reader_t *r)
{
	int i;
	for (i = 0; i < r->n_ctg; ++i) free(r->ctg);
	free(r->ctg);
	free(r->s.s);
	ks_destroy(r->ks);
	gzclose(r->fp);
	free(r);
}

int lt_cnt_read(lt_reader_t *r, lt_frag_t *f)
{
	int i, ret, dret;
	char *p, *q;
	if ((ret = ks_getuntil(r->ks, KS_SEP_LINE, &r->s, &dret)) < 0) return ret;
	for (p = q = r->s.s, i = 0;; ++p) {
		if (*p == 0 || *p == '\t') {
			*p = 0;
			if (i == 0) { // contig name
				if (r->n_ctg == 0 || strcmp(r->ctg[r->n_ctg-1], q) != 0) {
					if (r->n_ctg == r->m_ctg) {
						r->m_ctg = r->m_ctg? r->m_ctg<<1 : 4;
						r->ctg = (char**)realloc(r->ctg, r->m_ctg * sizeof(char*));
					}
					r->ctg[r->n_ctg++] = strdup(q);
				}
				f->ctg = r->n_ctg - 1;
			} else if (i == 1) {
				f->st = atoi(q);
			} else if (i == 2) {
				f->en = atoi(q);
			} else if (i == 3) {
				char *t;
				f->n_seg = strtol(q, &t, 10);
				f->n_frag = strtol(t + 1, &t, 10);
				assert(t[1] == '<' || t[1] == '|');
				f->lo = t[1] == '<'? 1 : 0;
				f->ro = t[2] == '>'? 1 : 0;
			}
			q = p + 1;
		}
	}
	return 0;
}

#include <unistd.h>

int main_count(int argc, char *argv[])
{
	int c;
	lt_reader_t *r;
	lt_frag_t f;
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti group <dedup.bam> | lianti count [options] -\n");
		return 1;
	}

	r = lt_cnt_open(argv[optind]);
	while (lt_cnt_read(r, &f) >= 0) {
	}
	lt_cnt_close(r);
	return 0;
}
