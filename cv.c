#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sam.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

#define NUM_BIN_SIZE 73

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);
uint64_t bed_totlen(void *_bed);

uint64_t bin_size_table[NUM_BIN_SIZE] = {
	1, 2, 3, 4, 5, 6, 7, 8, 9,
	10, 20, 30, 40, 50, 60, 70, 80, 90,
	100, 200, 300, 400, 500, 600, 700, 800, 900,
	1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
	10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
	100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
	1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000,
	10000000, 20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000,
	100000000
};

// update results after reading a new base pair, whose depth is n_plp
static void lt_cv_auto(uint64_t n_plp, uint64_t* _n_bp, double* current_value, double* n, double* sum, double* sqsum) {
	int i;
	++*_n_bp;
	for (i = 0; i < NUM_BIN_SIZE; ++i) {
		current_value[i] += n_plp;
		if (*_n_bp % bin_size_table[i] == 0) { // after a whole bin is read, store the 0th, 1st, and 2nd moments
			++n[i];
			sum[i] += current_value[i];
			sqsum[i] += current_value[i] * current_value[i];
			current_value[i] = 0;
		}
	}
}

int main_cv(int argc, char *argv[])
{
	int bed_input = 0;
	bam_plp_t plp;
	BGZF *fp;
	bam_hdr_t *h;
	const bam_pileup1_t *p;
	int i, c, n_plp, tid, pos, last_pos = -1, last_tid = -1;
	uint64_t bed_len = 0, n_bp = 0;
	double *current_value, *n, *sum, *sqsum;
	void *bed = 0;
	gzFile c_fp;
	kstream_t *ks;
	int dret, st = -1, en = -1;
	kstring_t *str_chr, *str_num;

	while ((c = getopt(argc, argv, "b:c")) >= 0) {
		if (c == 'b') {
			bed = bed_read(optarg);
			bed_len = bed_totlen(bed);
			fprintf(stderr, "[M::%s] total length in BED: %ld\n", __func__, (long)bed_len);
		} else if (c == 'c') {
			bed_input = 1;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti cv [-b bed] <aln.bam>\n");
		fprintf(stderr, "   or: lianti cv -c [-b bed] <count.bed.gz>\n");
		return 1;
	}

	current_value = (double*)calloc(NUM_BIN_SIZE, sizeof(double));
	n = (double*)calloc(NUM_BIN_SIZE, sizeof(double));
	sum = (double*)calloc(NUM_BIN_SIZE, sizeof(double));
	sqsum = (double*)calloc(NUM_BIN_SIZE, sizeof(double));
		
	if (bed_input) { // if input is a BED file ("-c") (for example, from "lianti count")
		// below are modified from bedidx.c and count.c
		c_fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
		str_chr = (kstring_t*)calloc(1, sizeof(kstring_t));
		str_num = (kstring_t*)calloc(1, sizeof(kstring_t));
		ks = ks_init(c_fp);
		while (ks_getuntil(ks, 0, str_chr, &dret) >= 0) { // read chr name
			for (i = 0; i < 3; ++i) {
				ks_getuntil(ks, 0, str_num, &dret);
				if (i == 0) st = atoi(str_num->s); // read region start
				else if (i == 1) en = atoi(str_num->s); // read region end
				else if (i == 2) n_plp = atoi(str_num->s); // read allele count
			}
			if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n'); // skip the rest of the line
			for (pos = st; pos < en; ++pos) { // check intersection with the mask bp-by-bp
                if (bed_overlap(bed, str_chr->s, pos, pos + 1) == 0)
                    continue;
				lt_cv_auto(n_plp, &n_bp, current_value, n, sum, sqsum);
			}
		}
		ks_destroy(ks);
		gzclose(c_fp);
		free(str_chr->s);
		free(str_num->s);
		free(str_chr);
		free(str_num);
		// above are modified from bedidx.c and count.c
	} else { // if input is a BAM file
		fp = bgzf_open(argv[optind], "r");
		h = bam_hdr_read(fp);

		plp = bam_plp_init((bam_plp_auto_f)bam_read1, fp);
		while ((p = bam_plp_auto(plp, &tid, &pos, &n_plp)) != 0) {
			// below are modified from develop branch of bam2depth.c (Jul 1, 2016)
	        while (tid > last_tid) {
	            if (last_tid >= 0) {
	                // Deal with remainder or entirety of last tid.
	                while (++last_pos < h->target_len[last_tid]) {
	                    // Horribly inefficient, but the bed API is an obfuscated black box.
	                    if (bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
	                        continue;
						lt_cv_auto(0, &n_bp, current_value, n, sum, sqsum);
	                }
	            }
	            last_tid++;
	            last_pos = -1;
	        }
	        // Deal with missing portion of current tid
	        while (++last_pos < pos) {
	            if (bed_overlap(bed, h->target_name[tid], last_pos, last_pos + 1) == 0)
	                continue;
				lt_cv_auto(0, &n_bp, current_value, n, sum, sqsum);
	        }
	        last_tid = tid;
	        last_pos = pos;	
			if (bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0)
				continue;
			lt_cv_auto(n_plp, &n_bp, current_value, n, sum, sqsum);
			// above are modified from develop branch of bam2depth.c (Jul 1, 2016)
		}
		bam_plp_destroy(plp);
		bam_hdr_destroy(h);
		bgzf_close(fp);
	}
	
	for (i = 0; i < NUM_BIN_SIZE; ++i) {
		printf("%ld\t%f\n", (long)bin_size_table[i], sqrt((sqsum[i] - sum[i] * sum[i] / n[i]) / n[i]) / sum[i] * n[i]); // calculate CV from the moments
	}
	
	free(current_value);
	free(n);
	free(sum);
	free(sqsum);
	if (bed) bed_destroy(bed);
	return 0;
}
