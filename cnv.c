#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include "kvec.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

#include "ksort.h"
KSORT_INIT_GENERIC(float)

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

msseg_t *mss_find_all(int n, const float *S, float min_sc, int *n_seg)
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

float mss_find_one(int n, const float *S)
{
	int i;
	float L, L_max;
	for (i = 0, L = L_max = 0.; i < n; ++i) {
		L += S[i];
		if (L < 0.) L = 0.;
		else if (L > L_max) L_max = L;
	}
	return L_max;
}

void mss_shuffle(int n, float *a)
{
	int i, j;
	for (i = n; i > 1; --i) {
		float tmp;
		j = (int)(drand48() * i);
		tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;
	}
}

/***********************
 * Gumbel distribution *
 ***********************/

void lt_gumbel_naive(int n, float *a, float x[2])
{
	int i;
	float mean, var, median;
	median = ks_ksmall_float(n, a, n/2);
	for (i = 0, mean = var = 0.; i < n; ++i) {
		mean += a[i];
		var += a[i] * a[i];
	}
	mean /= n;
	var /= n;
	x[1] = (mean - median) / (0.5772156649 - 0.3665129206);
	if (x[1] < 1.) x[1] = sqrt(6*var) / 3.1415926536;
	x[0] = mean - x[1] * 0.5772156649;
}

float lt_gumbel_update(int n, const float *a, float x[2])
{
	double s0, s1, s2, z[2], eps;
	int i;
	for (i = 0, s0 = s1 = s2 = 0; i < n; ++i) {
		double t = exp(-a[i] / x[1]);
		s0 += a[i];
		s1 += t;
		s2 += a[i] * t;
	}
	s0 /= n;
	z[0] = -x[1] * log(s1 / n);
	z[1] = s0 - s2 / s1;
	x[0] = fabs(x[0] - z[0]);
	x[1] = fabs(x[1] - z[1]);
	eps = x[0] > x[1]? x[0] : x[1];
	x[0] = z[0], x[1] = z[1];
	return eps;
}

double lt_gumbel_cdf(const float x[2], float z)
{
	return exp(-exp(-(z-x[0])/x[1]));
}

double lt_gumbel_ccdf(const float x[2], float z)
{
	double y;
	y = exp(-(z-x[0])/x[1]);
	return y > .001? 1. - lt_gumbel_cdf(x, z) : y * (1. - y * (.5 - y/6.));
}

double lt_gumbel_quantile(const float x[2], float p)
{
	return x[0] - x[1] * log(-log(p));
}

void lt_gumbel_est(int n, float *S, int n_perm, float x[2], int max_itr, float eps)
{
	int k;
	float *ev;
	ev = (float*)malloc(n * sizeof(float));
	for (k = 0; k < n_perm; ++k) {
		mss_shuffle(n, S);
		ev[k] = mss_find_one(n, S);
	}
	lt_gumbel_naive(n_perm, ev, x);
	for (k = 0; k < max_itr; ++k)
		if (lt_gumbel_update(n_perm, ev, x) < eps)
			break;
	if (isnan(x[1]) || isinf(x[1]) || k == max_itr) {
		fprintf(stderr, "[E::%s] failed to fit a Gumbel distribution. Abort!\n", __func__);
		abort();
	}
	free(ev);
}

/******************
 * CNV parameters *
 ******************/

typedef struct {
	int ploidy, n_perm, gumbel_max_itr;
	float pen_miss, pen_coef, gumbel_eps;
	float rep_thres;
} lt_cnvopt_t;

void lt_cnvopt_init(lt_cnvopt_t *opt)
{
	opt->ploidy = 2;
	opt->n_perm = 200;
	opt->pen_coef = 4.;
	opt->pen_miss = .001;
	opt->gumbel_max_itr = 1000;
	opt->gumbel_eps = .1;
	opt->rep_thres = .01;
}

/****************
 * Depth reader *
 ****************/

typedef struct {
	uint32_t e;
	uint8_t d[3], flt;
} lt_dp1_t;

typedef kvec_t(lt_dp1_t) lt_depth1_v;

typedef struct {
	char *name;
	lt_depth1_v d;
} lt_rawdp_t;

void lt_dp_destroy(int n, lt_rawdp_t *d)
{
	int i;
	for (i = 0; i < n; ++i) {
		free(d[i].name);
		free(d[i].d.a);
	}
	free(d);
}

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

lt_rawdp_t *lt_dp_read(const char *fn, int *n, const char *fn_gap)
{
	int dret;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	kvec_t(lt_rawdp_t) a = {0,0,0};
	lt_rawdp_t *r;
	void *gap = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	if (fn_gap) gap = bed_read(fn_gap);

	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int i, c, st;
		char *q, *p;
		lt_dp1_t t;
		for (i = 0, q = p = str.s;; ++p) {
			if (*p == 0 || *p == '\t') {
				c = *p, *p = 0;
				if (i == 0) {
					if (a.n == 0 || strcmp(q, a.a[a.n-1].name) != 0) {
						kv_pushp(lt_rawdp_t, a, &r);
						r->name = strdup(q);
						kv_init(r->d);
					}
				} else if (i == 1) {
					st = atoi(q);
				} else if (i == 2) {
					t.e = atoi(q);
				} else if (i >= 3 && i <= 5) {
					int x;
					x = atoi(q);
					t.d[i-3] = x < 255? x : 255;
				}
				if (c == 0) break;
				q = p + 1, ++i;
			}
		}
		if (i >= 5) {
			t.flt = gap && bed_overlap(gap, a.a[a.n-1].name, st, t.e)? 1 : 0;
			kv_push(lt_dp1_t, a.a[a.n-1].d, t);
		}
	}

	free(str.s);
	if (gap) bed_destroy(gap);
	ks_destroy(ks);
	gzclose(fp);
	*n = a.n;
	return a.a;
}

/*****************************
 * Estimate super-parameters *
 *****************************/

typedef struct {
	float pen_gain, gumbel_gain[2];
} lt_cnvpar_t;

static void gen_S(int n, const lt_dp1_t *d, int ploidy, float pen_nosig, float pen_miss, float *S)
{
	int i, l;
	for (i = l = 0; i < n; ++i) {
		const lt_dp1_t *p = &d[i];
		int len = p->e - (i? (p-1)->e : 0);
		if (p->flt) S[l++] = -pen_miss * len;
		else if (p->d[1] > ploidy) S[l++] = len;
		else if (p->d[0] < ploidy) S[l++] = -pen_nosig * len;
		else S[l++] = -pen_miss * len;
	}
}

void lt_dp_par(const lt_cnvopt_t *opt, int n_chr, const lt_rawdp_t *d, lt_cnvpar_t *par)
{
	int i, k, l, l_sig, l_nosig, tot;
	float *S;

	for (k = tot = 0; k < n_chr; ++k) tot += d[k].d.n;
	S = (float*)malloc(tot * sizeof(float));

	// gain
	l_sig = l_nosig = 0;
	for (k = 0; k < n_chr; ++k) {
		const lt_rawdp_t *dk = &d[k];
		for (i = 0; i < dk->d.n; ++i) {
			lt_dp1_t *p = &dk->d.a[i];
			int len = p->e - (i? (p-1)->e : 0);
			if (p->flt) continue;
			if (p->d[1] > opt->ploidy) l_sig += len;
			else if (p->d[0] < opt->ploidy) l_nosig += len;
		}
	}
	par->pen_gain = opt->pen_coef * l_sig / l_nosig;
	printf("GS\t%d\t%d\t%.3f\n", l_sig, l_nosig, par->pen_gain);

	for (k = l = 0; k < n_chr; ++k) {
		gen_S(d[k].d.n, d[k].d.a, opt->ploidy, par->pen_gain, opt->pen_miss, &S[l]);
		l += d[k].d.n;
	}
	lt_gumbel_est(tot, S, opt->n_perm, par->gumbel_gain, opt->gumbel_max_itr, opt->gumbel_eps);

	// loss

	// free
	free(S);
}

void lt_dp_cnv(const lt_cnvopt_t *opt, const lt_cnvpar_t *par, int n_chr, const lt_rawdp_t *dp)
{
	int max_len, k;
	float *S;
	for (k = max_len = 0; k < n_chr; ++k)
		max_len = max_len > dp[k].d.n? max_len : dp[k].d.n;
	S = (float*)malloc(max_len * sizeof(float));
	// gain
	for (k = 0; k < n_chr; ++k) {
		msseg_t *seg;
		int i, n_seg, n = dp[k].d.n;
		lt_dp1_t *d = dp[k].d.a;
		gen_S(n, d, opt->ploidy, par->pen_gain, opt->pen_miss, S);
		seg = mss_find_all(n, S, lt_gumbel_quantile(par->gumbel_gain, 1. - opt->rep_thres), &n_seg);
		for (i = 0; i < n_seg; ++i) {
			msseg_t *si = &seg[i];
			int j, en = d[si->en-1].e, st = si->st? d[si->st-1].e : 0;
			printf("GG\t%s\t%d\t%d\t%.2f\t%.3g\n", dp[k].name, st, en, si->sc, lt_gumbel_ccdf(par->gumbel_gain, si->sc));
			for (j = si->st; j < si->en; ++j)
				printf("GE\t%d\t%d\t%d\t%d\n", d[j].e, d[j].d[0], d[j].d[1], d[j].d[2]);
		}
		free(seg);
	}
	free(S);
}

int main_cnv(int argc, char *argv[])
{
	int c, n_chr;
	lt_cnvopt_t opt;
	lt_cnvpar_t par;
	lt_rawdp_t *dp;
	char *fn_gap = 0;

	lt_cnvopt_init(&opt);
	while ((c = getopt(argc, argv, "g:c:P:")) >= 0) {
		if (c == 'g') fn_gap = optarg;
		else if (c == 'c') opt.pen_coef = atof(optarg);
		else if (c == 'P') opt.rep_thres = atof(optarg);
		else if (c == 'p') opt.ploidy = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti cnv [options] <depth.bed>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -p INT      expected ploidy [%d]\n", opt.ploidy);
		fprintf(stderr, "  -P FLOAT    P-value threshold [%g]\n", opt.rep_thres);
		fprintf(stderr, "  -g FILE     places of assembly gaps in BED [null]\n");
		fprintf(stderr, "  -c FLOAT    penalty coefficient [%g]\n", opt.pen_coef);
		return 1;
	}

	dp = lt_dp_read(argv[optind], &n_chr, fn_gap);
	lt_dp_par(&opt, n_chr, dp, &par);
	lt_dp_cnv(&opt, &par, n_chr, dp);
	lt_dp_destroy(n_chr, dp);

	return 0;
}
