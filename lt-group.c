#include <stdlib.h>
#include <string.h>
#include "sam.h"
#include "kdq.h"

typedef struct {
	int l_ovlp;
	int max_seg;
} lt_opt_t;

typedef struct {
	uint32_t tid, st, en;
	uint32_t is_rev:1, is_closed:1;
	int n_frag;
} lt_group_t;

KDQ_INIT(lt_group_t)

typedef struct {
	kdq_t(lt_group_t) *q;
} lt_groups_t;

void lt_opt_init(lt_opt_t *opt)
{
	opt->l_ovlp = 9;
	opt->max_seg = 10000;
}

lt_groups_t *lt_grp_init(void)
{
	lt_groups_t *g;
	g = (lt_groups_t*)calloc(1, sizeof(lt_groups_t));
	g->q = kdq_init(lt_group_t);
	return g;
}

void lt_grp_destroy(lt_groups_t *g)
{
	kdq_destroy(lt_group_t, g->q);
	free(g);
}

void lt_grp_push(const lt_opt_t *opt, lt_groups_t *g, const bam_hdr_t *h, const bam1_t *b)
{
	int i, st, en, is_rev, is_closed = 0;
	const bam1_core_t *c = &b->core;
	const uint8_t *SA, *YD;

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
	YD = bam_aux_get(b, "YD");
	is_closed = YD? bam_aux2i(YD) : 0;
	is_closed = 0;

	while (kdq_size(g->q)) {
		lt_group_t *p = &kdq_first(g->q);
		if (p->tid != c->tid || p->en <= st) {
			printf("%s\t%d\t%d\t%d:%d\t%c\n", h->target_name[p->tid], p->st, p->en, p->n_frag, p->is_closed, "+-"[p->is_rev]);
			kdq_shift(lt_group_t, g->q);
		} else break;
	}
	for (i = 0; i < kdq_size(g->q); ++i) {
		lt_group_t *p = &kdq_at(g->q, i);
		int added = 0;
		if (p->is_closed) {
			if (is_closed) {
				if (st == p->st && en == p->en) added = 1;
			} else if (!is_rev) {
				if (st == p->st && en <= p->en) added = 1;
			} else {
				if (en == p->en) added = 1;
			}
		} else {
			if (is_closed) {
				if (!p->is_rev) {
					if (st == p->st && en >= p->en) added = 1, p->is_closed = 1;
				} else {
					if (en == p->en) added = 1, p->is_closed = 1;
				}
			} else {
				if (!p->is_rev) {
					if (st == p->st) added = 1, p->en = p->en > en? p->en : en;
				} else {
					if (en == p->en) added = 1;
				}
				if (!added) {
				}
			}
		}
		if (added) {
			++p->n_frag;
			break;
		}
	}
	if (i == kdq_size(g->q)) {
		lt_group_t *p;
		p = kdq_pushp(lt_group_t, g->q);
		p->tid = c->tid, p->st = st, p->en = en, p->is_rev = is_rev, p->is_closed = is_closed, p->n_frag = 1;
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
		lt_grp_push(&opt, g, h, b);
	bam_destroy1(b);
	lt_grp_destroy(g);

	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
