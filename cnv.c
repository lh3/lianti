#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "kvec.h"
#include "kseq.h"
//KSTREAM_DECLARE(gzFile, gzread)
KSTREAM_INIT(gzFile, gzread, 0x10000)

/*************************************
 * Find all maximal scoring segments *
 *************************************/

typedef struct {
	int st, en;
	float L, R;
	int pre;
} msseg_aux_t;

typedef struct {
	int st, en;
	float sc;
} msseg_t;

typedef kvec_t(msseg_t) msseg_v;
typedef kvec_t(msseg_aux_t) msseg_aux_v;

static void add_segs(msseg_v *ret, msseg_aux_v *seg, int min_sc)
{
	int i;
	for (i = 0; i < seg->n; ++i) {
		msseg_aux_t *p = &seg->a[i];
		if (p->R - p->L >= min_sc) {
			msseg_t *q;
			kv_pushp(msseg_t, *ret, &q);
			q->st = p->st, q->en = p->en, q->sc = p->R - p->L;
		}
	}
	seg->n = 0;
}

msseg_t *mssegs(int n, const float *S, int min_sc, int *n_seg)
{
	int i, j;
	float L;
	msseg_v ret = {0,0,0};
	msseg_aux_v seg = {0,0,0};
	msseg_aux_t t;

	for (i = L = 0; i < n;) {
		if (S[i] > 0) {
			int k;
			float R = L + S[i];
			for (k = i + 1; k < n && S[k] > 0.; ++k)
				R += S[k];
			t.st = i, t.en = k, t.L = L, t.R = R;
			while (1) {
				msseg_aux_t *p;
				for (j = seg.n - 1; j >= 0;) {
					p = &seg.a[j];
					if (p->L < t.L) break;
					j = p->pre >= 0? p->pre : j - 1;
				}
				if (j >= 0 && seg.a[j].R < t.R) {
					p = &seg.a[j];
					t.st = p->st, t.L = p->L, t.pre = p->pre;
					seg.n = j;
				} else {
					if (j < 0) add_segs(&ret, &seg, min_sc);
					t.pre = j;
					kv_push(msseg_aux_t, seg, t);
					break;
				}
			}
			L = R, i = k;
		} else L += S[i++];
	}
	add_segs(&ret, &seg, min_sc);
	free(seg.a);
	ret.a = (msseg_t*)realloc(ret.a, ret.n * sizeof(msseg_t));
	*n_seg = ret.n;
	return ret.a;
}

/******************/

typedef struct {
	int min_len, ploidy;
	float  ipen, dpen, min_sc;
} lt_cnvopt_t;

void lt_cnvopt_init(lt_cnvopt_t *opt)
{
	opt->ipen = 1.;
	opt->dpen = 1.;
	opt->min_sc = 1000.;
	opt->min_len = 10000;
	opt->ploidy = 2;
}

typedef struct {
	int e, d[3];
} lt_depth_t;

static void process(const lt_cnvopt_t *opt, const char *chr, int n, lt_depth_t *d)
{
	int i, n_seg;
	float *S;
	msseg_t *seg;
	S = (float*)malloc(n * sizeof(float));
	// gain
	for (i = 0; i < n; ++i) {
		lt_depth_t *di = &d[i];
		int st = i? d[i-1].e : 0;
		if (d[i].d[1] > opt->ploidy) S[i] = di->e - st;
		else if (d[i].d[0] <= opt->ploidy) S[i] = -opt->ipen * (di->e - st);
		else S[i] = 0.;
	}
	seg = mssegs(n, S, opt->min_sc, &n_seg);
	for (i = 0; i < n_seg; ++i) {
		msseg_t *si = &seg[i];
		int j, en = d[si->en-1].e, st = si->st? d[si->st-1].e : 0;
		if (en - st >= opt->min_len) {
			printf("G\t%s\t%d\t%d\t%.2f\n", chr, st, en, si->sc);
			for (j = si->st; j < si->en; ++j)
				printf("E\t%d\t%d\t%d\t%d\n", d[j].e, d[j].d[0], d[j].d[1], d[j].d[2]);
		}
	}
	// loss
	free(S);
}

int main_cnv(int argc, char *argv[])
{
	int c, dret;
	char *last_chr = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	kvec_t(lt_depth_t) dp = {0,0,0};
	lt_cnvopt_t opt;

	lt_cnvopt_init(&opt);
	while ((c = getopt(argc, argv, "s:l:i:d:")) >= 0) {
		if (c == 's') opt.min_sc = atof(optarg);
		else if (c == 'l') opt.min_len = atoi(optarg);
		else if (c == 'i') opt.ipen = atof(optarg);
		else if (c == 'd') opt.dpen = atof(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti cnv <depth.bed>\n");
		return 1;
	}

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);

	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int i = 0;
		char *q, *p;
		lt_depth_t t;
		for (q = p = str.s;; ++p) {
			if (*p == 0 || *p == '\t') {
				c = *p, *p = 0;
				if (i == 0) {
					if (last_chr && strcmp(last_chr, q) != 0) {
						process(&opt, last_chr, dp.n, dp.a);
						free(last_chr);
						dp.n = 0;
						last_chr = strdup(q);
					} else if (last_chr == 0) last_chr = strdup(q);
				} else if (i == 2) {
					t.e = atoi(q);
				} else if (i >= 3 && i <= 5) {
					t.d[i-3] = atoi(q);
				}
				if (c == 0) break;
				q = p + 1, ++i;
			}
		}
		if (i >= 5) kv_push(lt_depth_t, dp, t);
	}
	process(&opt, last_chr, dp.n, dp.a);
	free(last_chr);
	free(str.s);

	ks_destroy(ks);
	gzclose(fp);
	return 0;
}
