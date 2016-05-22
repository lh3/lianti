#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
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

/*************************
 * Brent 1D root finding *
 *************************/

#define BR_ITMAX 100
#define BR_EPS 3.0e-8
#define BR_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double brent_root(double (*func)(double, void*), double x1, double x2, double tol, int *err, void *data)
{
	int iter;
	double a = x1, b = x2, c = x2, d, e, min1, min2;
	double fa = (*func)(a, data), fb = (*func)(b, data), fc, p, q, r, s, tol1, xm;

	*err = 0;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		*err = 1;
		return 0.;
	}
	fc = fb, e = d = b - a;
	for (iter = 1; iter <= BR_ITMAX; ++iter) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
			c = a, fc = fa, e = d = b - a;
		if (fabs(fc) < fabs(fb)) {
			a = b, b = c, c = a;
			fa = fb, fb = fc, fc = fa;
		}
		tol1 = 2.0 * BR_EPS * fabs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if (a == c) {
				p = 2.0 * xm * s;
				q = 1.0 - s;
			} else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0 * xm * q - fabs(tol1 * q);
			min2 = fabs(e * q);
			if (2.0 * p < (min1 < min2 ? min1 : min2)) {
				e = d, d = p / q;
			} else d = xm, e = d;
		} else d = xm, e = d;
		a = b, fa = fb;
		if (fabs(d) > tol1) b += d;
		else b += BR_SIGN(tol1, xm);
		fb = (*func)(b, data);
	}
	*err = 2;
	return 0.;
}

/***********************
 * Gumbel distribution *
 ***********************/

typedef struct {
	int n;
	float *a;
	double mu;
} brent_aux_t;

static double gumbel_beta(double x, void *data)
{
	brent_aux_t *d = (brent_aux_t*)data;
	int i;
	double s0, s1, s2, beta;
	for (i = 0, s0 = s1 = s2 = 0.; i < d->n; ++i) {
		double t = exp(-d->a[i] / x);
		s0 += d->a[i];
		s1 += t;
		s2 += d->a[i] * t;
	}
	s0 /= d->n;
	beta = s0 - s2 / s1;
	d->mu = -x * log(s1 / d->n);
	return beta - x;
}

void lt_gumbel_est(int n, float *S, int n_perm, float x[2])
{
	brent_aux_t aux;
	int k, err;
	float *ev, x0 = 1000., x1, x2;
	double t;

	ev = (float*)malloc(n * sizeof(float));
	for (k = 0; k < n_perm; ++k) {
		mss_shuffle(n, S);
		ev[k] = mss_find_one(n, S);
	}
	aux.n = n_perm, aux.a = ev;
	while (1) {
		t = gumbel_beta(x0, &aux);
		if (isnan(t)) x0 *= 2.;
		else break;
	}
	if (t > 0.) {
		x1 = x0, x2 = x0 * 2.;
		while (1) {
			t = gumbel_beta(x2, &aux);
			if (t > 0.) x2 *= 2.;
			else break;
		}
	} else {
		x2 = x0, x1 = x0 / 2.;
		while (1) {
			t = gumbel_beta(x1, &aux);
			if (isnan(t)) x1 *= 1.414;
			else if (t < 0.) x1 /= 2;
			else break;
		}
	}
	x[1] = brent_root(gumbel_beta, x1, x2, 1e-3, &err, &aux);
	assert(err == 0);
	x[0] = aux.mu;
	free(ev);
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

/******************
 * CNV parameters *
 ******************/

typedef struct {
	int ploidy, n_perm, show_evidence, split_len;
	float pen_miss, pen_coef;
	float rep_thres;
} lt_cnvopt_t;

void lt_cnvopt_init(lt_cnvopt_t *opt)
{
	memset(opt, 0, sizeof(lt_cnvopt_t));
	opt->ploidy = 2;
	opt->n_perm = 200;
	opt->pen_coef = 4.;
	opt->pen_miss = .1;
	opt->rep_thres = 1e-4;
	opt->split_len = 1000;
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

lt_rawdp_t *lt_dp_read(const char *fn, int *n, const char *fn_gap, int split_len)
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
		int i, c, st = -1;
		char *q, *p;
		lt_dp1_t t;
		t.e = -1;
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
			if (t.flt == 0) {
				int len = t.e - st;
				lt_dp1_t t2 = t;
				while (len > split_len + (split_len>>1)) {
					t2.e = st + split_len;
					kv_push(lt_dp1_t, a.a[a.n-1].d, t2);
					len -= split_len;
					st += split_len;
				}
				t2.e = st + len;
				kv_push(lt_dp1_t, a.a[a.n-1].d, t2);
			} else kv_push(lt_dp1_t, a.a[a.n-1].d, t);
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

#define LT_GAIN 0
#define LT_LOSS1 1
#define LT_LOSS2 2

typedef struct {
	float penalty[3], gumbel[3][2];
} lt_cnvpar_t;

static inline int classify_signal(int type, const lt_dp1_t *p, int ploidy)
{
	int s = 0;
	if (p->flt) return 0;
	if (type == LT_GAIN) {
		if (p->d[1] > ploidy) s = 1;
		else if (p->d[0] <= ploidy) s = -1;
	} else if (type == LT_LOSS1) {
		if (p->d[2] < ploidy) s = 1;
		else if (p->d[2] >= ploidy) s = -1;
		else if (p->d[0] >= ploidy) s = -2;
	} else if (type == LT_LOSS2) {
		if (p->d[2] == 0) s = 1;
		else if (p->d[0] > 0) s = -4;
		else if (p->d[2] > 0) s = -8;
	}
	return s;
}

static void gen_S(const lt_cnvopt_t *opt, const lt_cnvpar_t *par, int type, int n, const lt_dp1_t *d, float *S)
{
	int i, l;
	float pen_nosig = par->penalty[type];
	float pen_miss = opt->pen_miss * par->penalty[type];
	for (i = l = 0; i < n; ++i) {
		const lt_dp1_t *p = &d[i];
		int len = p->e - (i? (p-1)->e : 0);
		int s = classify_signal(type, p, opt->ploidy);
		if (s == 1) S[l++] = len;
		else if (s < 0) S[l++] = s * pen_nosig * len;
		else S[l++] = -pen_miss * len;
	}
}

void lt_cnv_par(const lt_cnvopt_t *opt, int n_chr, const lt_rawdp_t *d, lt_cnvpar_t *par)
{
	int i, k, l, tot, type;
	float *S;

	for (k = tot = 0; k < n_chr; ++k) tot += d[k].d.n;
	S = (float*)malloc(tot * sizeof(float));

	for (type = 0; type < 3; ++type) {
		int64_t l_sig, l_nosig;
		l_sig = l_nosig = 0;
		for (k = 0; k < n_chr; ++k) {
			const lt_rawdp_t *dk = &d[k];
			for (i = 0; i < dk->d.n; ++i) {
				lt_dp1_t *p = &dk->d.a[i];
				int len = p->e - (i? (p-1)->e : 0);
				int s = classify_signal(type, p, opt->ploidy);
				if (s == 1) l_sig += len;
				else if (s < 0) l_nosig += len;
			}
		}
		par->penalty[type] = opt->pen_coef * l_sig / l_nosig;
		printf("%cS\t%ld\t%ld\t%.3f\n", "GLA"[type], (long)l_sig, (long)l_nosig, par->penalty[type]);
		for (k = l = 0; k < n_chr; ++k) {
			gen_S(opt, par, type, d[k].d.n, d[k].d.a, &S[l]);
			l += d[k].d.n;
		}
		lt_gumbel_est(tot, S, opt->n_perm, par->gumbel[type]);
		printf("%cP\t%.3f\t%.3f\t%.3f\n", "GLA"[type], par->gumbel[type][0], par->gumbel[type][1], lt_gumbel_quantile(par->gumbel[type], 1. - opt->rep_thres));
	}
	free(S);
}

void lt_cnv_call(const lt_cnvopt_t *opt, const lt_cnvpar_t *par, int n_chr, const lt_rawdp_t *dp)
{
	int max_len, k, type;
	float *S;
	for (k = max_len = 0; k < n_chr; ++k)
		max_len = max_len > dp[k].d.n? max_len : dp[k].d.n;
	S = (float*)malloc(max_len * sizeof(float));
	for (type = 0; type < 3; ++type) {
		for (k = 0; k < n_chr; ++k) {
			msseg_t *seg;
			int i, n_seg, n = dp[k].d.n;
			lt_dp1_t *d = dp[k].d.a;
			gen_S(opt, par, type, n, d, S);
			seg = mss_find_all(n, S, lt_gumbel_quantile(par->gumbel[type], 1. - opt->rep_thres), &n_seg);
			for (i = 0; i < n_seg; ++i) {
				msseg_t *si = &seg[i];
				int j, en = d[si->en-1].e, st = si->st? d[si->st-1].e : 0;
				printf("%c%c\t%s\t%d\t%d\t%.2f\t%.3g\n", "GLA"[type], "GLA"[type], dp[k].name, st, en, si->sc, lt_gumbel_ccdf(par->gumbel[type], si->sc));
				if (opt->show_evidence)
					for (j = si->st; j < si->en; ++j)
						printf("%cE\t%d\t%d\t%d\t%d\n", "GLA"[type], d[j].e, d[j].d[0], d[j].d[1], d[j].d[2]);
			}
			free(seg);
		}
	}
	free(S);
}

int main_cnv(int argc, char *argv[])
{
	int c, n_chr;
	lt_cnvopt_t opt;
	lt_cnvpar_t par;
	lt_rawdp_t *dp;

	lt_cnvopt_init(&opt);
	while ((c = getopt(argc, argv, "c:p:P:s:e")) >= 0) {
		if (c == 'P') opt.rep_thres = atof(optarg);
		else if (c == 'p') opt.ploidy = atoi(optarg);
		else if (c == 'e') opt.show_evidence = 1;
		else if (c == 'c') opt.pen_coef = atof(optarg);
		else if (c == 's') opt.split_len = atoi(optarg);
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: lianti cnv [options] <depth.bed> <assembly-gap.bed>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -p INT     expected ploidy [%d]\n", opt.ploidy);
		fprintf(stderr, "  -P FLOAT   P-value threshold [%g]\n", opt.rep_thres);
		fprintf(stderr, "  -c FLOAT   penalty coefficient [%g]\n", opt.pen_coef);
		fprintf(stderr, "  -s INT     split an interval if longer than INT [%d]\n", opt.split_len);
		return 1;
	}

	dp = lt_dp_read(argv[optind], &n_chr, argv[optind+1], opt.split_len);
	lt_cnv_par(&opt, n_chr, dp, &par);
	lt_cnv_call(&opt, &par, n_chr, dp);
	lt_dp_destroy(n_chr, dp);

	return 0;
}
