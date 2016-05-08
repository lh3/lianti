#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sam.h"
#include "kdq.h"
#include "kvec.h"
#include "ksort.h"
KSORT_INIT_GENERIC(int)

typedef struct {
	int l_ovlp;
	int max_seg;
	int min_frag;
	int fuzz, fuzz_ovlp;
	int no_merge;
} lt_opt_t;

typedef struct {
	int tid, st, en;
	int n_frag, n_seg;
	uint32_t is_rev:1, is_closed:1;

	int n, m;
	int *a;
} lt_group_t;

#define group_lt(a, b) ((a).st < (b).st)
KSORT_INIT(grp, lt_group_t, group_lt)

KDQ_INIT(lt_group_t)

typedef struct {
	kdq_t(lt_group_t) *q;
	int r_tid, r_max_en;

	int n, m;
	lt_group_t *a;
} lt_groups_t;

void lt_opt_init(lt_opt_t *opt)
{
	opt->l_ovlp = 9;
	opt->max_seg = 10000;
	opt->fuzz = 2;
	opt->fuzz_ovlp = 2;
	opt->min_frag = 2;
}

lt_groups_t *lt_grp_init(void)
{
	lt_groups_t *g;
	g = (lt_groups_t*)calloc(1, sizeof(lt_groups_t));
	g->q = kdq_init(lt_group_t);
	g->r_tid = -1;
	return g;
}

void lt_grp_destroy(lt_groups_t *g)
{
	kdq_destroy(lt_group_t, g->q);
	free(g->a);
	free(g);
}

void lt_grp_push_region(const lt_opt_t *opt, const bam_hdr_t *h, lt_groups_t *g, const lt_group_t *r)
{
	if (g->n && (r == 0 || r->tid != g->r_tid || r->st >= g->r_max_en)) {
		int i, j;
		ks_introsort(grp, g->n, g->a);
		// pass 1: forward-reverse merge
		for (i = 0; i < g->n; ++i) {
			lt_group_t *q = &g->a[i];
			if (q->n_frag == 0 || q->is_closed) continue;
			for (j = i + 1; j < g->n; ++j) {
				lt_group_t *p = &g->a[j];
				if (q->en <= p->st) break;
				if (p->n_frag == 0) continue;
				if (!q->is_rev && p->is_rev && q->en - p->en <= opt->fuzz && p->en - q->en <= opt->fuzz) { // merge
					q->n_frag += p->n_frag;
					q->is_closed = 1;
					q->en = q->en > p->en? q->en : p->en;
					p->n_frag = 0;
				} else if (q->is_rev != p->is_rev && p->st - q->st <= opt->fuzz) {
					if (q->en >= p->en) {
						if (q->is_rev) {
							q->n_frag += p->n_frag;
							q->is_closed = 1;
							p->n_frag = 0;
						}
					} else {
						if (p->is_rev) {
							p->n_frag += q->n_frag;
							p->is_closed = 1;
							q->n_frag = 0;
							break;
						}
					}
				}
			}
		}
		if (opt->no_merge) goto print_reg;
		/*/ pass 2: aggressive merge
		for (i = 0; i < g->n; ++i) {
			lt_group_t *q = &g->a[i];
			if (q->n_frag == 0 || q->is_rev || q->is_closed) continue;
			for (j = i + 1; j < g->n; ++j) {
				lt_group_t *p = &g->a[j];
				if (q->en <= p->st - opt->l_ovlp - 1) break;
				if (p->n_frag == 0 || !p->is_rev || p->en < q->en || p->is_closed) continue;
				q->n_frag += p->n_frag;
				q->en = p->en;
				++q->n_seg;
				q->is_closed = 1;
				p->n_frag = 0;
				break;
			}
		}*/
		// pass n: merge 9bp overlaps
		for (i = 0; i < g->n; ++i) {
			lt_group_t *q = &g->a[i];
			if (q->n_frag == 0) continue;
			for (j = i + 1; j < g->n; ++j) {
				lt_group_t *p = &g->a[j];
				if (q->en <= p->st) break;
				if (p->n_frag == 0) continue;
				if (q->en - p->st >= opt->l_ovlp - opt->fuzz_ovlp && q->en - p->st <= opt->l_ovlp + opt->fuzz_ovlp) { // TODO: better strategy: first count possible merges; don't merge if multiple
					q->n_frag += p->n_frag;
					q->en = p->en;
					++q->n_seg;
					p->n_frag = 0;
				}
			}
		}
print_reg:
		// print out
		for (i = 0; i < g->n; ++i) {
			lt_group_t *p = &g->a[i];
			if (p->n_frag)
				printf("%s\t%d\t%d\t%d:%d:%d:%d:%c\t%c\t%d\n", h->target_name[p->tid], p->st, p->en,
						p->st, p->en - p->st, p->n_frag, p->n_seg, p->is_closed? '*' : "+-"[p->is_rev], "+-"[p->is_rev], p->n_frag);
		}
		g->n = 0, g->r_tid = -1, g->r_max_en = 0;
	}
	if (r) {
		lt_group_t *t;
		kv_pushp(lt_group_t, *g, &t);
		memcpy(t, r, sizeof(lt_group_t));
		t->n_seg = 1;
		g->r_tid = t->tid;
		g->r_max_en = g->r_max_en > t->en? g->r_max_en : t->en;
	}
}

int lt_grp_segflt(int min_frag, lt_group_t *p)
{
	int flt;
	assert(p->n == p->n_frag);
	if (p->n >= min_frag) {
		int i, l1, l2 = 0, n2 = 0, s2;
		ks_introsort(int, p->n, p->a);
		l1 = p->a[p->n - min_frag];
		for (i = 1, s2 = 0, l2 = p->a[0]; i <= p->n; ++i) {
			if (i == p->n || p->a[i] != p->a[i-1]) {
				if (i - s2 > n2)
					n2 = i - s2, l2 = p->a[i-1];
				s2 = i;
			}
		}
		l1 = l1 > l2? l1 : l2;
		if (!p->is_rev) p->en = p->st + l1;
		else p->st = p->en - l1;
		flt = 0;
	} else flt = 1;
	free(p->a);
	p->a = 0, p->m = p->n = 0;
	return flt;
}

void lt_grp_push_read(const lt_opt_t *opt, lt_groups_t *g, const bam_hdr_t *h, const bam1_t *b)
{
	int i, st, en, is_rev;
	const bam1_core_t *c = &b->core;
	const uint8_t *SA;

	// compute st, en and rev
	if (c->flag & (BAM_FUNMAP|BAM_FDUP|BAM_FSUPP)) return;
	SA = bam_aux_get(b, "SA");
	if (SA) return; // TODO: we can relax this, in future
	if (c->flag & 2) {
		if (c->tid != c->mtid || c->isize > opt->max_seg || c->isize < 0)
			return;
		st = c->pos;
		en = st + c->isize;
		is_rev = (c->flag & BAM_FREAD1)? 0 : 1;
	} else {
		st = c->pos;
		en = st + bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
		is_rev = (c->flag & 16)? 1 : 0;
	}

	while (kdq_size(g->q)) {
		lt_group_t *p = &kdq_first(g->q);
		if (p->tid != c->tid || p->en <= st) {
			if (!lt_grp_segflt(opt->min_frag, p))
				lt_grp_push_region(opt, h, g, p);
			kdq_shift(lt_group_t, g->q);
		} else break;
	}
	for (i = 0; i < kdq_size(g->q); ++i) {
		lt_group_t *p = &kdq_at(g->q, i);
		int added = 0;
		if (p->is_rev != is_rev) continue;
		if (!p->is_rev) {
			if (st == p->st)
				added = 1, ++p->n_frag, p->en = en > p->en? en : p->en;
		} else {
			if (en == p->en)
				added = 1, ++p->n_frag;
		}
		if (added) {
			kv_push(int, *p, en - st);
			break;
		}
	}
	if (i == kdq_size(g->q)) {
		lt_group_t *p;
		p = kdq_pushp(lt_group_t, g->q);
		p->tid = c->tid, p->st = st, p->en = en, p->is_rev = is_rev, p->n_frag = 1;
		p->n = 0, p->m = 4;
		p->a = (int*)malloc(p->m * sizeof(int));
		p->a[p->n++] = en - st;
	}
}

#include <unistd.h>

int main(int argc, char *argv[])
{
	int c;
	lt_opt_t opt;
	BGZF *fp;
	bam_hdr_t *h;
	bam1_t *b;
	lt_groups_t *g;

	lt_opt_init(&opt);
	while ((c = getopt(argc, argv, "l:n:f:M")) >= 0) {
		if (c == 'l') opt.l_ovlp = atoi(optarg);
		else if (c == 'n') opt.min_frag = atoi(optarg);
		else if (c == 'f') opt.fuzz = atoi(optarg);
		else if (c == 'M') opt.no_merge = 1;
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lt-group [options] <in.bam>\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bgzf_open(argv[optind], "r") : bgzf_dopen(fileno(stdin), "r");
	h = bam_hdr_read(fp);

	g = lt_grp_init();
	b = bam_init1();
	while (bam_read1(fp, b) >= 0)
		lt_grp_push_read(&opt, g, h, b);
	lt_grp_push_region(&opt, h, g, 0);
	bam_destroy1(b);
	lt_grp_destroy(g);

	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
