#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include "kvec.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

#define MAX_HIST 10
static uint64_t cnt_hist[MAX_HIST+1];

extern void ks_introsort_int(size_t, int[]);

typedef struct {
	int ctg;
	uint32_t lo:1, st:31;
	uint32_t ro:1, en:31;
	int mq;
	int n_seg, n_frag;
} lt_frag_t;

#include "kdq.h"
KDQ_INIT(lt_frag_t)

typedef struct {
	int no_merge;
	int min_frag, min_frag2;
	int min_mq;
} lt_copt_t;

void lt_copt_init(lt_copt_t *opt)
{
	memset(opt, 0, sizeof(lt_copt_t));
	opt->min_frag = 5;
	opt->min_frag2 = 10;
	opt->min_mq = 40;
}

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
	for (i = 0; i < r->n_ctg; ++i) free(r->ctg[i]);
	free(r->ctg);
	free(r->s.s);
	ks_destroy(r->ks);
	gzclose(r->fp);
	free(r);
}

int lt_cnt_read(lt_reader_t *r, lt_frag_t *f)
{
	int i, ret, c, dret;
	char *p, *q;
	if ((ret = ks_getuntil(r->ks, KS_SEP_LINE, &r->s, &dret)) < 0) return ret;
	for (p = q = r->s.s, i = 0;; ++p) {
		if (*p == 0 || *p == '\t') {
			c = *p, *p = 0;
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
			} else if (i == 6) {
				f->mq = atoi(q);
			}
			if (c == 0) break;
			++i, q = p + 1;
		}
	}
	return 0;
}

typedef struct {
	kdq_t(lt_frag_t) *q;
	int last_ctg, last_pos;
} lt_cntbuf_t;

lt_cntbuf_t *lt_buf_init(void)
{
	lt_cntbuf_t *b;
	b = (lt_cntbuf_t*)calloc(1, sizeof(lt_cntbuf_t));
	b->q = kdq_init(lt_frag_t);
	return b;
}

void lt_buf_destroy(lt_cntbuf_t *b)
{
	kdq_destroy(lt_frag_t, b->q);
	free(b);
}

static void clear_up_to(const lt_copt_t *opt, lt_cntbuf_t *b, int end, char *const* ctg)
{
	int s = b->last_pos;
	kdq_t(lt_frag_t) *q = b->q;
	while (kdq_size(q) && s < end) {
		int i, s2 = end, d = 0, d2 = 0, d3 = 0;
		for (i = 0; i < kdq_size(q); ++i) {
			lt_frag_t *f = &kdq_at(q, i);
			if (s >= f->st && s < f->en) { // overlapping the counting position
				if (f->mq >= opt->min_mq) ++d;
				if (f->mq >= opt->min_mq && f->n_frag >= opt->min_frag2) ++d2;
				++d3;
				if (f->en <= end)
					s2 = s2 < f->en? s2 : f->en;
			} else if (s < f->st) { // start after the counting position
				s2 = s2 < f->st? s2 : f->st;
			}
		}
		if (s2 != INT_MAX) {
			printf("%s\t%d\t%d\t%d\t%d\t%d\n", ctg[b->last_ctg], s, s2, d, d2, d3);
			cnt_hist[d < MAX_HIST? d : MAX_HIST] += s2 - s;
		}
		s = s2;
		while (kdq_size(q) && kdq_first(q).en <= s)
			kdq_shift(lt_frag_t, q);
	}
	if (end != INT_MAX) {
		cnt_hist[0] += end - s;
		if (s < end) printf("%s\t%d\t%d\t0\t0\t0\n", ctg[b->last_ctg], s, end);
	}
	b->last_pos = end;
}

static int test_merge(lt_cntbuf_t *b, lt_frag_t *f)
{
	if (f->ctg == b->last_ctg && f->lo) {
		int i, max = 0, max_i = -1;
		kdq_t(lt_frag_t) *q = b->q;
		for (i = 0; i < kdq_size(q); ++i) {
			lt_frag_t *g = &kdq_at(q, i);
			if (g->ro && f->st < g->en && f->en >= g->en) {
				if (g->en - f->st > max) max = g->en - f->st, max_i = i;
			}
		}
		if (max > 0) {
			lt_frag_t *g = &kdq_at(q, max_i);
			assert(f->en >= g->en && f->st < g->en && f->st >= g->st);
			g->ro = f->ro, g->en = f->en;
			++g->n_seg, g->n_frag += f->n_frag;
			return 1;
		}
	}
	return 0;
}

void lt_buf_push(const lt_copt_t *opt, lt_cntbuf_t *b, lt_frag_t *f, char *const* ctg)
{
	if (f) {
		lt_frag_t *p;
		if (!opt->no_merge && test_merge(b, f)) return;
		if (f->n_frag < opt->min_frag) return;
		if (f->ctg != b->last_ctg) {
			clear_up_to(opt, b, INT_MAX, ctg);
			b->last_ctg = f->ctg;
			b->last_pos = 0;
		}
		clear_up_to(opt, b, f->st, ctg);
		b->last_pos = f->st;
		p = kdq_pushp(lt_frag_t, b->q);
		memcpy(p, f, sizeof(lt_frag_t));
	} else clear_up_to(opt, b, INT_MAX, ctg);
}

#include <unistd.h>

int main_count(int argc, char *argv[])
{
	int c, i;
	lt_copt_t opt;
	lt_reader_t *r;
	lt_frag_t f;
	lt_cntbuf_t *b;

	lt_copt_init(&opt);
	while ((c = getopt(argc, argv, "Mn:q:")) >= 0) {
		if (c == 'M') opt.no_merge = 1;
		else if (c == 'q') opt.min_mq = atoi(optarg);
		else if (c == 'n') {
			char *q;
			opt.min_frag = strtol(optarg, &q, 10);
			opt.min_frag2 = *q == ','? atoi(q+1) : opt.min_frag*2;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti group <dedup.bam> | lianti count [options] -\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -n INT1[,INT2]  ignore fragments consisting of <INT reads [%d,INT1*2]\n", opt.min_frag);
		fprintf(stderr, "  -M              do not merge open-ended overlapping fragments\n");
		fprintf(stderr, "  -q INT          min RMS mapping quality [%d]\n\n", opt.min_mq);
		fprintf(stderr, "Output:\n");
		fprintf(stderr, "  chr  start  end  depth{mapQ>=%d&&nReads>=%d} depth{mapQ>=%d&&nReads>=%d} depthAll\n", opt.min_mq, opt.min_frag, opt.min_mq, opt.min_frag2);
		return 1;
	}

	r = lt_cnt_open(argv[optind]);
	b = lt_buf_init();
	while (lt_cnt_read(r, &f) >= 0)
		lt_buf_push(&opt, b, &f, r->ctg);
	lt_buf_push(&opt, b, 0, r->ctg);
	lt_buf_destroy(b);
	lt_cnt_close(r);
	for (i = 0; i <= MAX_HIST; ++i)
		fprintf(stderr, "H\t%d\t%ld\n", i, (long)cnt_hist[i]);
	return 0;
}
