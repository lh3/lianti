#include <stdlib.h>
#include <string.h>
#include "sam.h"
#include "kdq.h"

typedef struct {
	int l_ovlp;
	int max_seg;
	int fuzz;
} lt_opt_t;

typedef struct {
	uint32_t tid:31, is_rev:1;
	int st, en;
	int n_frag, n_seg;
} lt_group_t;

KDQ_INIT(lt_group_t)

typedef struct {
	kdq_t(lt_group_t) *q;
	kdq_t(lt_group_t) *r;
	int r_tid, r_max_en;
} lt_groups_t;

void lt_opt_init(lt_opt_t *opt)
{
	opt->l_ovlp = 9;
	opt->max_seg = 10000;
	opt->fuzz = 2;
}

lt_groups_t *lt_grp_init(void)
{
	lt_groups_t *g;
	g = (lt_groups_t*)calloc(1, sizeof(lt_groups_t));
	g->q = kdq_init(lt_group_t);
	g->r = kdq_init(lt_group_t);
	g->r_tid = -1;
	return g;
}

void lt_grp_destroy(lt_groups_t *g)
{
	kdq_destroy(lt_group_t, g->r);
	kdq_destroy(lt_group_t, g->q);
	free(g);
}

void lt_grp_push_region(const lt_opt_t *opt, const bam_hdr_t *h, lt_groups_t *g, const lt_group_t *r)
{
	if (kdq_size(g->r) && (r == 0 || r->tid != g->r_tid || r->st >= g->r_max_en)) {
		int i, j;
		// pass 1: merge obvious forward-reverse fragments
		for (i = 0; i < kdq_size(g->r); ++i) {
			lt_group_t *q = &kdq_at(g->r, i);
			for (j = i + 1; j < kdq_size(g->r); ++j) {
				lt_group_t *p = &kdq_at(g->r, j);
				if (q->en <= p->st) break;
				if (p->n_frag == 0) continue;
				if (!q->is_rev && p->is_rev && q->en == p->en) { // merge
					q->n_frag += p->n_frag;
					p->n_frag = 0;
				} else if (q->is_rev != p->is_rev && q->st == p->st) {
					if (q->en > p->en) {
						if (q->is_rev) {
							q->n_frag += p->n_frag;
							p->n_frag = 0;
						}
					} else {
						if (p->is_rev) {
							p->n_frag += q->n_frag;
							q->n_frag = 0;
							break;
						}
					}
				}
			}
		}
		// pass n: merge 9bp overlaps
		for (i = 0; i < kdq_size(g->r); ++i) {
			lt_group_t *q = &kdq_at(g->r, i);
			if (q->n_frag == 0) continue;
			for (j = i + 1; j < kdq_size(g->r); ++j) {
				lt_group_t *p = &kdq_at(g->r, j);
				if (q->en <= p->st) break;
				if (p->n_frag == 0) continue;
				if (q->en - p->st == opt->l_ovlp || q->en - p->st == opt->l_ovlp + 1 || q->en - p->st == opt->l_ovlp - 1) { // TODO: better strategy: first count possible merges; don't merge if multiple
					q->n_frag += p->n_frag;
					q->en = p->en;
					++q->n_seg;
					p->n_frag = 0;
					break;
				}
			}
		}
		// fuzzy merge
		for (i = 0; i < kdq_size(g->r); ++i) {
			lt_group_t *q = &kdq_at(g->r, i);
			if (q->n_frag == 0) continue;
			for (j = i + 1; j < kdq_size(g->r); ++j) {
				lt_group_t *p = &kdq_at(g->r, j);
				if (q->en <= p->st) break;
				if (p->n_frag == 0) continue;
				if (q->n_seg == p->n_seg && ((q->st - p->st <= opt->fuzz && p->st - q->st <= opt->fuzz) || (q->en - p->en <= opt->fuzz && p->en - q->en <= opt->fuzz))) {
					q->n_frag += p->n_frag;
					q->en = q->en > p->en? q->en : p->en;
					p->n_frag = 0;
					break;
				}
			}
		}
		// print out
		while (kdq_size(g->r)) {
			lt_group_t *p = kdq_shift(lt_group_t, g->r);
			if (p->n_frag)
				printf("%s\t%d\t%d\t%d:%d\t%c\t%d\n", h->target_name[p->tid], p->st, p->en, p->n_frag, p->n_seg+1, "+-"[p->is_rev], p->n_frag);
		}
		g->r_max_en = 0;
	}
	kdq_push(lt_group_t, g->r, *r);
	kdq_last(g->r).n_seg = 1;
	g->r_tid = r->tid;
	g->r_max_en = g->r_max_en > r->en? g->r_max_en : r->en;
}

void lt_grp_push_read(const lt_opt_t *opt, lt_groups_t *g, const bam_hdr_t *h, const bam1_t *b)
{
	int i, st, en, is_rev;
	const bam1_core_t *c = &b->core;
	const uint8_t *SA;

	// compute st, en and rev
	if (c->flag & 4) return; // unmapped
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
			//printf("%s\t%d\t%d\t*\t%c\t%d\n", h->target_name[p->tid], p->st, p->en, "+-"[p->is_rev], p->n_frag);
			lt_grp_push_region(opt, h, g, p);
			kdq_shift(lt_group_t, g->q);
		} else break;
	}
	for (i = 0; i < kdq_size(g->q); ++i) {
		lt_group_t *p = &kdq_at(g->q, i);
		int added = 0;
		if (!p->is_rev) {
			if (st == p->st) added = 1, ++p->n_frag, p->en = p->en > en? p->en : en;
		} else {
			if (en == p->en) added = 1, ++p->n_frag;
		}
		if (added) break;
	}
	if (i == kdq_size(g->q)) {
		lt_group_t *p;
		p = kdq_pushp(lt_group_t, g->q);
		p->tid = c->tid, p->st = st, p->en = en, p->is_rev = is_rev, p->n_frag = 1;
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
	while ((c = getopt(argc, argv, "l:")) >= 0) {
		if (c == 'l') opt.l_ovlp = atoi(optarg);
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
	bam_destroy1(b);
	lt_grp_destroy(g);

	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}