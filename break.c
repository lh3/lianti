#include <stdio.h>
#include <stdlib.h>
#include "sam.h"

typedef struct {
	const char *seq;
	int pos, mapq, nm;
	int qs, qe, rs, re, rev;
} side_t;

int main_break(int argc, char *argv[])
{
	BGZF *fp;
	bam_hdr_t *h;
	bam1_t *b;
	int c, min_mapq = 20;

	while ((c = getopt(argc, argv, "q:")) >= 0) {
		if (c == 'q') min_mapq = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti break [options] <in.bam>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -q INT   min mapping quality [%d]\n", min_mapq);
		return 1;
	}

	fp = bgzf_open(argv[optind], "r");
	h = bam_hdr_read(fp);

	b = bam_init1();
	while (bam_read1(fp, b) >= 0) {
		const bam1_core_t *c = &b->core;
		const uint8_t *SA = 0;
		char *sa, *p;
		int i, n_semicolon = 0;
		side_t s[2], t;
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
				s[i].seq = h->target_name[c->tid];
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
				s[i].seq = sa;
				for (p = sa; *p != ','; ++p);
				*p = 0;
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
			if (s[i].mapq < min_mapq) break;
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
		printf("%s\t%d\t%c\t%d\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%d\t%d\n", s[0].seq, s[0].re, "+-"[s[0].rev], s[0].mapq, s[0].qe - s[0].qs, s[0].nm,
				s[1].seq, s[1].rs, "+-"[s[1].rev], s[1].mapq, s[1].qe - s[1].qs, s[1].nm, s[1].qs - s[0].qe);
	}
	bam_destroy1(b);

	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
