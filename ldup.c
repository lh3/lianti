#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "sam.h"

typedef struct {
	bam1_t *b;
	uint8_t is_pe, is_left, full_ovlp;
} elem_t;

#include "kdq.h"
KDQ_INIT(elem_t)

#include "khash.h"
KHASH_SET_INIT_INT64(64)
KHASH_SET_INIT_STR(s)

#define AUX_REALLOC_SIZE 1024

static uint64_t lt_n_frags_noBC, lt_n_dups_noBC, lt_n_frags_BC, lt_n_dups_BC;

static inline uint64_t X31_hash_string(const char *s)
{
	uint64_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

static khash_t(64) *process(kdq_t(elem_t) *q, BGZF *fp, khash_t(s) *marked, khash_t(64) *aux)
{
	int i, absent;
	if (kh_n_buckets(aux) < kdq_size(q) * 4 && kh_n_buckets(aux) <= AUX_REALLOC_SIZE) {
		kh_clear(64, aux);
	} else {
		kh_destroy(64, aux);
		aux = kh_init(64);
	}
	for (i = 0; i < kdq_size(q); ++i) {
		elem_t *e = &kdq_at(q, i);
		bam1_t *b = e->b;
		char *qname = bam_get_qname(b);
		if (b->core.flag&0x800) continue;
		if ((b->core.flag&1) && (b->core.flag&(BAM_FREAD1|BAM_FREAD2))) { // PE
			e->is_pe = 1;
			e->full_ovlp = (b->core.tid == b->core.mtid && b->core.pos == b->core.mpos);
			e->is_left = 0;
			if (b->core.tid < b->core.mtid) e->is_left = 1;
			else if (b->core.tid == b->core.mtid) {
				if (b->core.pos < b->core.mpos) e->is_left = 1;
				else if (b->core.pos == b->core.mpos && (b->core.flag&BAM_FREAD1)) e->is_left = 1;
			}
		} else { // SE
			e->is_left = 1;
			e->is_pe = e->full_ovlp = 0;
		}
		if (e->is_left) {
			const uint8_t *BC = 0;
			uint64_t key;
			khint_t k;
			BC = bam_aux_get(b, "BC");
			if (BC) ++lt_n_frags_BC;
			else ++lt_n_frags_noBC;
			key = BC? X31_hash_string(bam_aux2Z(BC)) : 0;
			k = kh_put(64, aux, key, &absent);
			if (!absent) {
				if (BC) ++lt_n_dups_BC;
				else ++lt_n_dups_noBC;
				b->core.flag |= BAM_FDUP;
				if (e->is_pe) kh_put(s, marked, strdup(qname), &absent);
			}
		}
	}
	for (i = 0; i < kdq_size(q); ++i) {
		elem_t *e = &kdq_at(q, i);
		bam1_t *b = e->b;
		char *qname = bam_get_qname(b);
		if (e->is_pe && !e->is_left) {
			khint_t k;
			k = kh_get(s, marked, qname);
			if (k != kh_end(marked)) {
				b->core.flag |= BAM_FDUP;
				free((char*)kh_key(marked, k));
				kh_del(s, marked, k);
			}
		}
	}
	while (kdq_size(q)) {
		elem_t *e;
		e = kdq_shift(elem_t, q);
		bam_write1(fp, e->b);
		bam_destroy1(e->b);
	}
	return aux;
}

int main_ldup(int argc, char *argv[])
{
	int c, clevel = -1, ret;
	int last_tid = -1, last_pos = -1;
	BGZF *fpr, *fpw;
	bam_hdr_t *h;
	bam1_t *b;
	khash_t(s) *marked;
	khash_t(64) *aux;
	kdq_t(elem_t) *q;
	khint_t k;

	while ((c = getopt(argc, argv, "l:")) >= 0) {
		if (c == 'l') clevel = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti ldup [options] <aln.bam>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    zlib compression level [zlib default]\n");
		return 1;
	}

	fpr = strcmp(argv[optind], "-")? bgzf_open(argv[optind], "r") : bgzf_dopen(fileno(stdin), "r");
	h = bam_hdr_read(fpr);
	if (clevel >= 0 && clevel <= 9) {
		char mode[3];
		sprintf(mode, "w%d", clevel);
		fpw = bgzf_dopen(fileno(stdout), mode);
	} else fpw = bgzf_dopen(fileno(stdout), "w");
	bgzf_mt(fpw, 3, 256);
	bam_hdr_write(fpw, h);
	aux = kh_init(64);
	marked = kh_init(s);
	q = kdq_init(elem_t);

	b = bam_init1();
	while ((ret = bam_read1(fpr, b)) >= 0) {
		elem_t *e;
		b->core.flag &= ~BAM_FDUP;
		if (b->core.tid != last_tid || b->core.pos != last_pos) {
			if (last_tid >= 0 && last_pos >= 0)
				aux = process(q, fpw, marked, aux);
			last_tid = b->core.tid, last_pos = b->core.pos;
		}
		if (b->core.tid < 0) break;
		e = kdq_pushp(elem_t, q);
		e->b = bam_init1();
		bam_copy1(e->b, b);
	}
	aux = process(q, fpw, marked, aux);
	if (ret >= 0) {
		do {
			bam_write1(fpw, b);
		} while (bam_read1(fpr, b) >= 0);
	}
	bam_destroy1(b);

	kdq_destroy(elem_t, q);
	fprintf(stderr, "[M::%s] %ld+%ld fragments; %ld+%ld duplicates; %d unpaired reads\n", __func__,
			(long)lt_n_frags_BC, (long)lt_n_frags_noBC, (long)lt_n_dups_BC, (long)lt_n_dups_noBC, kh_size(marked));
	for (k = 0; k < kh_end(marked); ++k) 
		if (kh_exist(marked, k)) free((char*)kh_key(marked, k));
	kh_destroy(s, marked);
	kh_destroy(64, aux);
	bgzf_close(fpw);
	bam_hdr_destroy(h);
	bgzf_close(fpr);
	return 0;
}
