#include <stdlib.h>
#include <stdio.h>
#include "sam.h"

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);
uint64_t bed_totlen(void *_bed);

#define MAX_DEPTH 10000

int main_lorenz(int argc, char *argv[])
{
	bam_plp_t plp;
	BGZF *fp;
	bam_hdr_t *h;
	const bam_pileup1_t *p;
	int i, c, n_plp, tid, pos, step = 1000;
	uint64_t *cnt, bed_len = 0, sum_partial = 0, cov = 0, tot = 0, tot_partial = 0;
	void *bed = 0;

	while ((c = getopt(argc, argv, "b:s:")) >= 0) {
		if (c == 'b') {
			bed = bed_read(optarg);
			bed_len = bed_totlen(bed);
			fprintf(stderr, "[M::%s] total length in BED: %ld\n", __func__, (long)bed_len);
		} else if (c == 's') step = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: lianti lorenz [-b bed] [-s step=%d] <aln.bam>\n", step);
		return 1;
	}

	cnt = (uint64_t*)calloc(MAX_DEPTH + 1, sizeof(uint64_t));
	fp = bgzf_open(argv[optind], "r");
	h = bam_hdr_read(fp);

	plp = bam_plp_init((bam_plp_auto_f)bam_read1, fp);
	while ((p = bam_plp_auto(plp, &tid, &pos, &n_plp)) != 0) {
		if (bed_overlap(bed, h->target_name[tid], pos, pos + 1))
			++cnt[n_plp < MAX_DEPTH? n_plp : MAX_DEPTH];
	}
	for (i = 1; i <= MAX_DEPTH; ++i) cov += cnt[i];
	cnt[0] = bed_len - cov;
	for (i = 0; i <= MAX_DEPTH; ++i) tot += cnt[i] * i;
	bam_plp_destroy(plp);

	printf("%.4f\t%.4f\n", 0., 0.);
	for (i = 0, sum_partial = tot_partial = 0; i <= MAX_DEPTH; ++i) {
		if (cnt[i] <= step) {
			sum_partial += cnt[i], tot_partial += i * cnt[i];
			printf("%.4f\t%.4f\n", (double)sum_partial / bed_len, (double)tot_partial / tot);
		} else {
			uint64_t rest = cnt[i];
			while (rest) {
				int x = rest < step? rest : step;
				sum_partial += x, tot_partial += i * x;
				printf("%.4f\t%.4f\n", (double)sum_partial / bed_len, (double)tot_partial / tot);
				rest -= x;
			}
		}
	}

	bam_hdr_destroy(h);
	bgzf_close(fp);
	free(cnt);
	if (bed) bed_destroy(bed);
	return 0;
}
