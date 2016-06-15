#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include "sam.h"
#include "kvec.h"
#include "ksort.h"

typedef struct {
	int tid, pos, mapq, nm;
	int qs, qe, rs, re, rev;
} side_t;

typedef struct {
	int64_t mid;
	uint64_t pos[2];
	int qgap;
	int dir:16, rdrev:16;
} break_t;

#define sv1_lt(a, b) ((a).mid < (b).mid)
KSORT_INIT(sv1, break_t, sv1_lt)

#define sv2_lt(a, b) ((a).pos[0] < (b).pos[0])
KSORT_INIT(sv2, break_t, sv2_lt)

int main_sv(int argc, char *argv[])
{
	BGZF *fp;
	bam_hdr_t *h;
	bam1_t *b;
	int c, i, min_mapq = 40, max_nm = 4, print_bp = 0, max_gap = 50, min_cnt = 3;
	uint64_t *off;
	kvec_t(break_t) a = {0,0,0};

	while ((c = getopt(argc, argv, "pq:m:g:n:")) >= 0) {
		if (c == 'q') min_mapq = atoi(optarg);
		else if (c == 'm') max_nm = atoi(optarg);
		else if (c == 'p') print_bp = 1;
		else if (c == 'g') max_gap = atoi(optarg);
		else if (c == 'n') min_cnt = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti sv [options] <in.bam>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -q INT   min mapping quality [%d]\n", min_mapq);
		fprintf(stderr, "  -m INT   max NM [%d]\n", max_nm);
		fprintf(stderr, "  -g INT   max gap [%d]\n", max_gap);
		fprintf(stderr, "  -n INT   min count [%d]\n", min_cnt);
		fprintf(stderr, "  -p       output break points, not the clustered SV calls\n");
		return 1;
	}

	fp = bgzf_open(argv[optind], "r");
	h = bam_hdr_read(fp);

	off = (uint64_t*)calloc(h->n_targets + 1, 8);
	for (i = 0; i < h->n_targets; ++i) off[i+1] = off[i] + h->target_len[i];

	b = bam_init1();
	while (bam_read1(fp, b) >= 0) {
		const bam1_core_t *c = &b->core;
		const uint8_t *SA = 0;
		char *sa, *p;
		int i, n_semicolon = 0;
		side_t s[2], t;
		int64_t mid;
		if ((c->flag & (BAM_FUNMAP|BAM_FSUPP|BAM_FQCFAIL|BAM_FSECONDARY|BAM_FDUP)) || c->tid < 0) continue;
		SA = bam_aux_get(b, "SA");
		if (SA == 0) continue;
		sa = bam_aux2Z(SA);
		for (p = sa; *p; ++p)
			if (*p == ';') ++n_semicolon;
		if (n_semicolon != 1) continue;

		for (i = 0; i < 2; ++i) {
			int k, clip[2], ql, rl;
			clip[0] = clip[1] = ql = rl = 0;
			if (i == 0) {
				const uint8_t *NM;
				const uint32_t *cigar;
				cigar = bam_get_cigar(b);
				s[i].tid = c->tid;
				s[i].rs = c->pos;
				s[i].rev = (c->flag&BAM_FREVERSE)? 1 : 0;
				for (k = 0; k < c->n_cigar; ++k) {
					int op = bam_cigar_op(cigar[k]);
					int len = bam_cigar_oplen(cigar[k]);
					if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) clip[k?1:0] = len;
					else if (op == BAM_CMATCH) ql += len, rl += len;
					else if (op == BAM_CINS) ql += len;
					else if (op == BAM_CDEL || op == BAM_CREF_SKIP) rl += len;
				}
				s[i].mapq = c->qual;
				s[i].nm = ((NM = bam_aux_get(b, "NM")) != 0)? bam_aux2i(NM) : -1;
			} else {
				int n_op = 0;
				for (p = sa; *p != ','; ++p);
				*p = 0;
				s[i].tid = bam_name2id(h, sa);
				assert(s[i].tid >= 0 && s[i].tid < h->n_targets);
				s[i].rs = strtol(p+1, &p, 10) - 1;
				s[i].rev = p[1] == '+'? 0 : 1;
				for (p += 3; *p && *p != ','; ++p, ++n_op) {
					int len;
					len = strtol(p, &p, 10);
					if (*p == 'H' || *p == 'S') clip[n_op?1:0] = len;
					else if (*p == 'M') ql += len, rl += len;
					else if (*p == 'I') ql += len;
					else if (*p == 'D' || *p == 'N') rl += len;
				}
				s[i].mapq = strtol(p+1, &p, 10);
				s[i].nm = strtol(p+1, &p, 10);
			}
			if (s[i].mapq < min_mapq || s[i].nm > max_nm) break;
			if (clip[0] == 0 && clip[1] == 0) break;
			if (s[i].rev) {
				s[i].qs = clip[1];
				s[i].re = s[i].rs;
				s[i].rs = s[i].re + rl;
			} else {
				s[i].qs = clip[0];
				s[i].re = s[i].rs + rl;
			}
			s[i].qe = s[i].qs + ql;
		}
		if (i != 2) continue;
		if (s[0].qs > s[1].qs) t = s[0], s[0] = s[1], s[1] = t;
		mid = ((off[s[0].tid] + s[0].re) + (off[s[1].tid] + s[1].rs)) >> 1;
		mid = (!s[0].rev && s[0].rs + off[s[0].tid] < s[1].rs + off[s[1].tid]) || (s[0].rev && s[0].rs + off[s[0].tid] > s[1].rs + off[s[1].tid])? mid : mid + off[h->n_targets];
		if (s[0].rev != s[1].rev) mid = -mid;
		if (!print_bp) {
			break_t *p;
			uint64_t tmp;
			kv_pushp(break_t, a, &p);
			p->mid = mid;
			p->pos[0] = (uint64_t)s[0].tid<<32 | s[0].re;
			p->pos[1] = (uint64_t)s[1].tid<<32 | s[1].rs;
			p->qgap = s[1].qs - s[0].qe;
			p->dir = (!!s[0].rev)<<1 | (!!s[1].rev);
			p->rdrev = 0;
			if (p->pos[0] > p->pos[1]) {
				tmp = p->pos[0], p->pos[0] = p->pos[1], p->pos[1] = tmp;
				p->dir = ((p->dir&1)^1)<<1 | (p->dir>>1^1);
				p->rdrev = 1;
			}
		} else printf("%s\t%d\t%c\t%d\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%d\t%d\t%lld\n", h->target_name[s[0].tid], s[0].re, "+-"[s[0].rev], s[0].mapq, s[0].qe - s[0].qs, s[0].nm,
					h->target_name[s[1].tid], s[1].rs, "+-"[s[1].rev], s[1].mapq, s[1].qe - s[1].qs, s[1].nm, s[1].qs - s[0].qe, (long long)mid);
	}
	bam_destroy1(b);

	if (!print_bp) {
		int start;
		ks_introsort(sv1, a.n, a.a);
		for (start = 0, i = 1; i <= a.n; ++i) {
			if (i == a.n || a.a[i].mid - a.a[i-1].mid > max_gap) {
				if (i - start >= min_cnt) {
					int j, subst;
					ks_introsort(sv2, i - start, &a.a[start]);
					for (subst = start, j = start + 1; j <= i; ++j) {
						if (j == i || a.a[j].pos[0]>>32 != a.a[j-1].pos[0]>>32 || a.a[j].pos[1]>>32 != a.a[j-1].pos[1]>>32 || a.a[j].pos - a.a[j-1].pos > max_gap) {
							int k, type;
							int64_t pos[2], qgap = 0;
							break_t *p = &a.a[subst];
							for (k = subst, pos[0] = pos[1] = 0; k < j; ++k) {
								pos[0] += (uint32_t)a.a[k].pos[0];
								pos[1] += (uint32_t)a.a[k].pos[1];
								qgap += a.a[k].qgap;
							}
							pos[0] = (int)((double)pos[0] / (j - subst) + .499);
							pos[1] = (int)((double)pos[1] / (j - subst) + .499);
							qgap = (int)((double)qgap / (j - subst) + .499);
							type = p->mid >= 0 && p->mid < off[h->n_targets] && p->pos[0]>>32 == p->pos[1]>>32? 'G' : p->pos[0]>>32 == p->pos[1]>>32? 'S' : 'T';
							printf("SV\t%s\t%d\t%c\t%s\t%d\t%c\t%c", h->target_name[p->pos[0]>>32], (uint32_t)pos[0], "+-"[p->dir>>1&1],
									h->target_name[p->pos[1]>>32], (uint32_t)pos[1], "+-"[p->dir&1], type);
							if (type == 'G') printf("\t%d", (int)(qgap - (pos[1] - pos[0])));
							else printf("\t.");
							printf("\t%d\n", j - subst);
							for (k = subst; k < j; ++k) {
								break_t *q = &a.a[k];
								printf("RD\t%s\t%d\t%c\t%s\t%d\t%c\t%d\t%c\n", h->target_name[q->pos[0]>>32], (uint32_t)q->pos[0], "+-"[q->dir>>1&1],
										h->target_name[q->pos[1]>>32], (uint32_t)q->pos[1], "+-"[q->dir&1], q->qgap, "+-"[q->rdrev]);
							}
							printf("//\n");
							subst = j;
						}
					}
				}
				start = i;
			}
		}
		free(a.a);
	}

	free(off);
	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
