#include <stdlib.h>
#include <string.h>
#include "sam.h"

typedef struct {
	int l_ovlp;
} lt_opt_t;

typedef struct {
	uint32_t s, e;
} lt_intv_t;

typedef struct {
	uint32_t st, en;
	int is5, n, m;
	lt_intv_t *a;
} lt_group_t;

typedef struct {
	uint32_t n, m;
	lt_group_t *a;
} lt_groups_t;

void lt_opt_init(lt_opt_t *opt)
{
	opt->l_ovlp = 9;
}

lt_groups_t *lt_grp_init(void)
{
	return (lt_groups_t*)calloc(1, sizeof(lt_groups_t));
}

void lt_grp_destroy(lt_groups_t *g)
{
	free(g->a);
	free(g);
}

void lt_grp_push(lt_groups_t *g, const bam_hdr_t *h, const bam1_t *b)
{
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
		lt_grp_push(g, h, b);
	bam_destroy1(b);

	lt_grp_destroy(g);
	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
