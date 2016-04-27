const char *lt_bind = "GGGAGATGTGTATAAGAGACAG"; // including the leading GGG
const char *lt_bind_rev = "CTGTCTCTTATACACATCTCCC";
const char *lt_promoter = "GAACAGAATTTAATACGACTCACTATA";

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include "kvec.h"

/******************
 * K-mer matching *
 ******************/

#include "khash.h"
KHASH_SET_INIT_INT64(s64)
typedef khash_t(s64) lt_seqcloud1_t;

typedef struct {
	int l;
	uint64_t s;
	lt_seqcloud1_t *mm;
	lt_seqcloud1_t *ins, *del; // generated but not used for now
} lt_seqcloud_t;

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static lt_seqcloud_t *lt_sc_init(void)
{
	lt_seqcloud_t *sc;
	sc = (lt_seqcloud_t*)calloc(1, sizeof(lt_seqcloud_t));
	sc->mm  = kh_init(s64);
	sc->ins = kh_init(s64);
	sc->del = kh_init(s64);
	return sc;
}

void lt_sc_destroy(lt_seqcloud_t *sc)
{
	kh_destroy(s64, sc->mm);
	kh_destroy(s64, sc->ins);
	kh_destroy(s64, sc->del);
	free(sc);
}

void lt_sc_add_core(lt_seqcloud_t *sc, uint64_t s)
{
	int i, absent;
	sc->s = s = s & ((1ULL<<sc->l*2) - 1);
	// mismatch
	for (i = 0; i < sc->l; ++i) {
		int i2 = i * 2, a, c = s>>i2&3;
		for (a = 1; a < 4; ++a) {
			uint64_t x = (s & ~(3ULL << i2)) | (uint64_t)((a+c)&3) << i2;
			kh_put(s64, sc->mm, x, &absent);
		}
	}

	// The following is actually not used for now

	// insertion
	for (i = 1; i < sc->l - 1; ++i) {
		int i2 = i * 2, a;
		uint64_t mask = (1ULL<<i2) - 1;
		for (a = 0; a < 4; ++a) {
			uint64_t x = (s & mask) | (uint64_t)a << i2 | s >> i2 << (i2 + 2);
			kh_put(s64, sc->ins, x, &absent);
		}
	}
	// deletion
	for (i = 1; i < sc->l - 1; ++i) {
		uint64_t x = (s & ((1ULL<<i*2) - 1)) | s >> (i+1)*2 << i*2;
		kh_put(s64, sc->del, x, &absent);
	}
}

lt_seqcloud_t *lt_sc_gen(const char *s)
{
	lt_seqcloud_t *sc;
	uint64_t x = 0;
	int i;
	sc = lt_sc_init();
	sc->l = strlen(s);
	for (i = 0; s[i] && i < sc->l; ++i) {
		int c = seq_nt4_table[(uint8_t)s[i]];
		if (c > 3) {
			lt_sc_destroy(sc);
			return 0;
		}
		x = x << 2 | c;
	}
	lt_sc_add_core(sc, x);
	return sc;
}

typedef struct {
	uint32_t pos:30, type:2;
} lt_sc_hit_t;

typedef struct { uint32_t n, m; lt_sc_hit_t *a; } lt_sc_hit_v;

lt_sc_hit_v lt_sc_test(const lt_seqcloud_t *sc, const char *seq)
{
	int i, l;
	uint64_t x = 0, mask = (1ULL << sc->l*2) - 1;
	lt_sc_hit_v hits = {0,0,0};
	lt_sc_hit_t tmp;
	for (i = l = 0; seq[i]; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x = (x << 2 | c) & mask;
			if (++l >= sc->l) {
				if (x == sc->s) {
					tmp.pos = i - (sc->l - 1), tmp.type = 0;
					kv_push(lt_sc_hit_t, hits, tmp);
				} else if (kh_get(s64, sc->mm, x) != kh_end(sc->mm)) {
					tmp.pos = i - (sc->l - 1), tmp.type = 1;
					kv_push(lt_sc_hit_t, hits, tmp);
				}
			}
		} else l = 0, x = 0;
	}
	return hits;
}

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: sc <kmer> <seq>\n");
		return 1;
	}
	lt_seqcloud_t *sc = lt_sc_gen(argv[1]);
	lt_sc_hit_v hits = lt_sc_test(sc, argv[2]);
	int i;
	for (i = 0; i < hits.n; ++i)
		printf("%d\t%d\n", hits.a[i].pos, hits.a[i].type);
	return 0;
}
