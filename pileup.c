// This piece of code is modified from samtools/bam2depth.c
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include "sam.h"
#include "faidx.h"
#include "ksort.h"

const char *hts_parse_reg(const char *s, int *beg, int *end);
void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

typedef struct {     // auxiliary data structure
	BGZF *fp;        // the file handler
	hts_itr_t *itr;  // NULL if a region not specified
	const bam_hdr_t *h;
	int min_mapQ, min_len; // mapQ filter; length filter
	int min_supp_len;
	float div_coef;
	void *bed;       // bedidx if not NULL
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->itr? bam_itr_next(aux->fp, aux->itr, b) : bam_read1(aux->fp, b);
	if (ret < 0) return ret;
	if (b->core.tid < 0) b->core.flag |= BAM_FUNMAP;
	if (!(b->core.flag&BAM_FUNMAP)) {
		if ((int)b->core.qual < aux->min_mapQ) {
			b->core.flag |= BAM_FUNMAP;
		} else if (aux->min_len > 0 || aux->min_supp_len > 0 || aux->bed) {
			int k, qlen = 0, tlen = 0;
			const char *chr = aux->h->target_name[b->core.tid];
			const uint32_t *cigar = bam_get_cigar(b);
			for (k = 0; k < b->core.n_cigar; ++k) { // compute the query length in the alignment
				int op = bam_cigar_op(cigar[k]);
				int oplen = bam_cigar_oplen(cigar[k]);
				if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP)
					qlen += oplen;
				if (bam_cigar_type(op)&2)
					tlen += oplen;
			}
			if (qlen < aux->min_len) b->core.flag |= BAM_FUNMAP;
			if (qlen < aux->min_supp_len && (b->core.flag&BAM_FSUPP)) b->core.flag |= BAM_FUNMAP;
			if (aux->bed && !(b->core.flag&BAM_FUNMAP) && !bed_overlap(aux->bed, chr, b->core.pos, b->core.pos + tlen))
				b->core.flag |= BAM_FUNMAP;
		}
		if (!(b->core.flag&BAM_FUNMAP) && aux->div_coef < 1.) {
			uint8_t *NM;
			int nm, k, n_gaps = 0, n_opens = 0, n_matches = 0;
			const uint32_t *cigar = bam_get_cigar(b);
			if ((NM = bam_aux_get(b, "NM")) == 0) return ret;
			if ((nm = bam_aux2i(NM)) == 0) return ret;
			for (k = 0; k < b->core.n_cigar; ++k) {
				int op = bam_cigar_op(cigar[k]);
				int l = bam_cigar_oplen(cigar[k]);
				if (op == BAM_CMATCH) n_matches += l;
				else if (op == BAM_CINS || op == BAM_CDEL) ++n_opens, n_gaps += l;
			}
			if (n_gaps <= nm) {
				int x = (nm - n_gaps) + n_opens, q;
				double expected = (n_matches + n_gaps) * aux->div_coef;
				double y = 1., p = 1.;
				if (x < expected) return ret;
				for (k = 1; k < x; ++k)
					y *= expected / k, p += y;
				p = 1. - p * exp(-expected);
				p = p > 1e-6? -4.343 * log(p) : 60.;
				q = (int)(p + .499);
				b->core.qual = b->core.qual > q? b->core.qual - q : 0;
				if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
			}
		}
	}
	return ret;
}

typedef struct {
	uint32_t is_skip:1, is_rev:1, b:4, q:8, k:18; // b=base, q=quality, k=allele id
	int indel; // <0: deleteion; >0: insertion
	uint32_t lt_pos;
	uint64_t hash;
	uint64_t pos; // i<<32|j: j-th read of the i-th sample
} allele_t;

#define allele_lt(a, b) ((a).hash < (b).hash || ((a).hash == (b).hash && (a).indel < (b).indel))
KSORT_INIT(allele, allele_t, allele_lt)

#define allelelt_lt(a, b) ((a).lt_pos < (b).lt_pos || ((a).lt_pos == (b).lt_pos && allele_lt(a, b)))
KSORT_INIT(allelelt, allele_t, allelelt_lt)

static inline allele_t pileup2allele(const bam_pileup1_t *p, int min_baseQ, uint64_t pos, int ref, int trim_len)
{ // collect allele information given a pileup1 record
	allele_t a;
	int i;
	const bam1_core_t *c = &p->b->core;
	const uint8_t *seq = bam_get_seq(p->b);
	a.k = (1<<18)-1; // this will be set in count_alleles()
	a.q = bam_get_qual(p->b)[p->qpos];
	a.is_rev = bam_is_rev(p->b);
	a.is_skip = (p->is_del || p->is_refskip || a.q < min_baseQ);
	if (p->qpos < trim_len || p->b->core.l_qseq - p->qpos < trim_len) a.is_skip = 1;
	a.indel = p->indel;
	a.b = a.hash = bam_seqi(seq, p->qpos);
	a.pos = pos;
	a.lt_pos = UINT32_MAX;
	if (bam_aux_get(p->b, "BC") != 0) {
		if (!(c->flag & BAM_FSUPP)) {
			uint32_t pos5, is_rev = c->flag&BAM_FREVERSE? 1 : 0;
			pos5 = is_rev? c->pos + bam_cigar2rlen(c->n_cigar, bam_get_cigar(p->b)) - 1 : c->pos;
			if (c->flag & BAM_FPAIRED) {
				if (c->flag & BAM_FPROPER_PAIR)
					a.lt_pos = (c->flag & BAM_FREAD1? pos5 : pos5 + c->isize) << 1 | is_rev;
			} else a.lt_pos = pos5<<1 | is_rev;
		}
	}
	if (p->indel > 0) // compute the hash for the insertion
		for (i = 0; i < p->indel; ++i)
			a.hash = (a.hash<<4) + a.hash + bam_seqi(seq, p->qpos + i + 1);
	a.hash = a.hash << 1 >> 1;
	if (p->indel != 0 || a.b != ref || ref == 15) // the highest bit tells whether it is a reference allele or not
		a.hash |= 1ULL<<63;
	return a;
}

static inline void print_allele(const bam_pileup1_t *p, int l_ref, const char *ref, int pos, int max_del, int is_vcf)
{ // print the allele. The format depends on is_vcf.
	const uint8_t *seq = bam_get_seq(p->b);
	int i, rest = max_del;
	putchar(seq_nt16_str[bam_seqi(seq, p->qpos)]);
	if (p->indel > 0) {
		if (!is_vcf) printf("+%d", p->indel);
		for (i = 1; i <= p->indel; ++i)
			putchar(seq_nt16_str[bam_seqi(seq, p->qpos + i)]);
	} else if (p->indel < 0) {
		if (!is_vcf) {
			printf("%d", p->indel);
			for (i = 1; i <= -p->indel; ++i)
				putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
		} else rest -= -p->indel, pos += -p->indel;
	}
	if (is_vcf)
		for (i = 1; i <= rest; ++i)
			putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
}

static void lt_drop_reads(int n, allele_t *a)
{
	int sti, i, j;
	ks_introsort(allelelt, n, a);
	for (sti = 0, i = 1; i <= n; ++i) {
		if (i == n || a[i-1].lt_pos != a[i].lt_pos) { // change of fragment
			int max_indel = 0, max = 0, max2 = 0, stj;
			uint64_t max_hash = 0;
			for (stj = sti, j = sti + 1; j <= i; ++j) {
				if (a[j].indel != a[j-1].indel || a[j].hash != a[j-1].hash) { // change of allele
					int cnt = j - sti;
					if (cnt > max) max2 = max, max = cnt, max_indel = a[j].indel, max_hash = a[j].hash;
					else if (cnt > max2) max2 = cnt;
					stj = j;
				}
			}
			if (max == max2) continue;
			for (j = sti; j <= i; ++j)
				if (a[j].indel != max_indel || a[j].hash != max_hash) // drop non-optimal reads
					a[j].is_skip = 1;
			sti = i;
			if (a[i].lt_pos == UINT32_MAX) break;
		}
	}
}

typedef struct {
	int n_a, n_alleles, max_del; // n_a: #reads used to compute quality sum; max_del: max deletion length
	int tot_dp, max_dp, n_cnt, max_cnt;
	allele_t *a; // allele of each read, of size $n_a
	int *cnt_strand, *cnt_supp; // cnt_strand: count of supporting reads on both strands; cnt_supp: sum of both strands
	int *support; // support across entire $a. It points to the last "row" of cnt_q.

	int len, max_len;
	char *seq;
	int *depth;
} paux_t;

static void count_alleles(paux_t *pa, int n, int qual_as_depth)
{
	allele_t *a = pa->a;
	int i, j;
	a[0].k = 0; // the first allele is given allele id 0
	pa->max_del = a[0].indel < 0? -a[0].indel : 0;
	for (i = pa->n_alleles = 1; i < pa->n_a; ++i) {
		if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) // change of allele
			++pa->n_alleles;
		a[i].k = pa->n_alleles - 1;
		pa->max_del = pa->max_del > -a[i].indel? pa->max_del : -a[i].indel; // max deletion
	}
	// collect per-BAM counts
	pa->n_cnt = pa->n_alleles * (n + 1);
	if (pa->n_cnt > pa->max_cnt) { // expand the arrays if necessary
		pa->max_cnt = pa->n_cnt;
		kroundup32(pa->max_cnt);
		pa->cnt_strand = (int*)realloc(pa->cnt_strand, pa->max_cnt * 2 * sizeof(int));
		pa->cnt_supp = (int*)realloc(pa->cnt_supp, pa->max_cnt * sizeof(int));
	}
	memset(pa->cnt_strand, 0, pa->n_cnt * 2 * sizeof(int));
	memset(pa->cnt_supp, 0, pa->n_cnt * sizeof(int));
	pa->support = pa->cnt_supp + pa->n_alleles * n; // points to the last row of cnt_q
	for (i = 0; i < pa->n_a; ++i) { // compute counts and sums of qualities
		int d = qual_as_depth? a[i].q : 1;
		j = (a[i].pos>>32)*pa->n_alleles + a[i].k;
		pa->cnt_strand[j<<1|a[i].is_rev] += d;
		pa->cnt_supp[j] += d;
		pa->support[a[i].k] += d;
	}
}

static void write_fa(paux_t *a, const char *name, int beg, float max_dev, int l_ref)
{
	int i, n_pos, max_dp;
	uint64_t sum_dp;
	double avg_dp, max_dp_real;
	if (l_ref == 0) l_ref = INT_MAX;
	for (i = 0, sum_dp = 0, n_pos = 0; i < a->len; ++i)
		if (a->seq[i] != 'n' && a->seq[i] != 'N')
			++n_pos, sum_dp += a->depth[i];
	avg_dp = (double)sum_dp/n_pos;
	max_dp_real = avg_dp + max_dev * sqrt(avg_dp);
	max_dp = max_dp_real > 0x7fffffff? 0x7fffffff : (int)(max_dp_real + .499); // to avoid integer overflow
	if (max_dp < 255) {
		for (i = 0; i < a->len; ++i)
			if (a->depth[i] > max_dp)
				a->seq[i] = tolower(a->seq[i]);
	}
	printf(">%s", name);
	if (beg > 0) printf(":%d", beg + 1);
	for (i = 0; i < a->len && i < l_ref; ++i) {
		if (i%60 == 0) putchar('\n');
		putchar(a->seq[i]);
	}
	if (l_ref < INT_MAX)
		for (; i < l_ref; ++i) {
			if (i%60 == 0) putchar('\n');
			putchar('N');
		}
	putchar('\n');
	fprintf(stderr, "[M::%s] average depth for contig '%s': %.2f\n", __func__, name, avg_dp);
}

int main_pileup(int argc, char *argv[])
{
	int i, j, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, l_ref = 0, min_support = 1, min_supp_len = 0, lt_mode = 0;
	int qual_as_depth = 0, is_vcf = 0, var_only = 0, show_2strand = 0, is_fa = 0, majority_fa = 0, rand_fa = 0, trim_len = 0, char_x = 'X';
	int last_tid, last_pos;
	float max_dev = 3.0, div_coef = 1.;
	const bam_pileup1_t **plp;
	char *ref = 0, *reg = 0, *chr_end; // specified region
	char *fname = 0; // reference fasta
	faidx_t *fai = 0;
	bam_hdr_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	paux_t aux;
	bam_mplp_t mplp;
	void *bed = 0;

	// parse the command line
	while ((n = getopt(argc, argv, "r:q:Q:l:f:dvcCS:Fs:D:V:uRMb:T:x:L")) >= 0) {
		if (n == 'f') { fname = optarg; fai = fai_load(fname); }
		else if (n == 'b') bed = bed_read(optarg);
		else if (n == 'l') min_len = atoi(optarg); // minimum query length
		else if (n == 'r') reg = strdup(optarg);   // parsing a region requires a BAM header
		else if (n == 'Q') baseQ = atoi(optarg);   // base quality threshold
		else if (n == 'q') mapQ = atoi(optarg);    // mapping quality threshold
		else if (n == 's') min_support = atoi(optarg);
		else if (n == 'd') qual_as_depth = 1;
		else if (n == 'S') min_supp_len = atoi(optarg);
		else if (n == 'v') var_only = 1;
		else if (n == 'V') div_coef = atof(optarg);
		else if (n == 'c') is_vcf = var_only = 1;
		else if (n == 'C') show_2strand = 1;
		else if (n == 'D') max_dev = atof(optarg), is_fa = 1;
		else if (n == 'F') is_fa = 1;
		else if (n == 'M') majority_fa = is_fa = 1;
		else if (n == 'R') rand_fa = is_fa = 1;
		else if (n == 'T') trim_len = atoi(optarg);
		else if (n == 'x') char_x = toupper(*optarg);
		else if (n == 'L') lt_mode = 1;
		else if (n == 'u') {
			baseQ = 3; mapQ = 20; qual_as_depth = 1;
			min_supp_len = 300; min_support = 5; div_coef = .01;
		}
	}
	if (min_support < 1) min_support = 1;
	if (is_fa && is_vcf) {
		fprintf(stderr, "[E::%s] option -F cannot be used with -c\n", __func__);
		return 1;
	}
	if (majority_fa && rand_fa) {
		fprintf(stderr, "[E::%s] option -M and -R can't be applied at the same time\n", __func__);
		return 1;
	}
	if (is_fa) var_only = 0;
	if (is_fa && !majority_fa && !rand_fa && min_support <= 1)
		fprintf(stderr, "[W::%s] with option -F, setting a reasonable -s is highly recommended.\n", __func__);
	if (is_vcf && fai == 0) {
		fprintf(stderr, "[E::%s] with option -c, the reference genome must be provided.\n", __func__);
		return 1;
	}
	if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   pileup [options] in1.bam [in2.bam [...]]\n\n");
        fprintf(stderr, "Options: -f FILE    reference genome [null]\n");
		fprintf(stderr, "         -r STR     region [null]\n");
		fprintf(stderr, "         -b FILE    BED or position list file to include [null]\n");
		fprintf(stderr, "         -Q INT     minimum base quality [%d]\n", baseQ);
		fprintf(stderr, "         -q INT     minimum mapping quality [%d]\n", mapQ);
		fprintf(stderr, "         -l INT     minimum query length [%d]\n", min_len);
		fprintf(stderr, "         -S INT     minimum supplementary alignment length [0]\n");
		fprintf(stderr, "         -V FLOAT   ignore queries with per-base divergence >FLOAT [1]\n");
		fprintf(stderr, "         -T INT     ignore bases within INT-bp from either end of a read [0]\n");
		fprintf(stderr, "         -d         base quality as depth\n");
		fprintf(stderr, "         -s INT     drop alleles with depth<INT [%d]\n", min_support);
		fprintf(stderr, "         -L         Lianti mode: choose the best allele among fragments with the same 5'-start\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "         -v         show variants only\n");
		fprintf(stderr, "         -c         output in the VCF format (force -v)\n");
		fprintf(stderr, "         -C         show count of each allele on both strands\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "         -F         output the consensus in FASTA\n");
		fprintf(stderr, "         -M         majority-allele FASTA (majfa; force -F)\n");
		fprintf(stderr, "         -R         random-allele FASTA (randfa; force -F)\n");
		fprintf(stderr, "         -x CHAR    character for bases identical to the reference [%c]\n", char_x);
		fprintf(stderr, "         -D FLOAT   soft mask if sumQ > avgSum+FLOAT*sqrt(avgSum) (force -F) [%.2f]\n", max_dev);
		fprintf(stderr, "\n");
		fprintf(stderr, "         -u         unitig calling mode (-d -V.01 -S300 -q20 -Q3 -s5)\n");
		fprintf(stderr, "\n");
		return 1;
	}

	// initialize the auxiliary data structures
	n = argc - optind; // the number of BAMs on the command line
	if (is_fa && n > 1) {
		fprintf(stderr, "[W::%s] with option -F, only the first input file is used.\n", __func__);
		n = 1;
	}
	srand48(11);
	data = (aux_t**)calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	if (reg) {
		chr_end = (char*)hts_parse_reg(reg, &beg, &end);
		ref = fai? fai_fetch(fai, reg, &l_ref) : 0;
	} else chr_end = 0;

	// load the index or put the file position at the right place
	last_tid = -1; last_pos = beg - 1;
	for (i = 0; i < n; ++i) {
		bam_hdr_t *htmp;
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = bgzf_open(argv[optind+i], "r"); // open BAM
		data[i]->min_mapQ = mapQ;                     // set the mapQ filter
		data[i]->min_len  = min_len;                  // set the qlen filter
		data[i]->div_coef = div_coef;
		data[i]->min_supp_len = min_supp_len;
		data[i]->bed = bed;
		htmp = bam_hdr_read(data[i]->fp);             // read the BAM header
		if (i == 0 && chr_end) {
			char c = *chr_end;
			*chr_end = 0;
			last_tid = tid = bam_name2id(htmp, reg);
			*chr_end = c;
		}
		if (i) bam_hdr_destroy(htmp); // if not the 1st BAM, trash the header
		else h = htmp; // keep the header of the 1st BAM
		if (tid >= 0) { // if a region is specified and parsed successfully
			hts_idx_t *idx = bam_index_load(argv[optind+i]); // load the index
			data[i]->itr = bam_itr_queryi(idx, tid, beg, end); // set the iterator
			hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
		}
		data[i]->h = h;
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = (int*)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = (const bam_pileup1_t**)calloc(n, sizeof(const bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
	memset(&aux, 0, sizeof(paux_t));
	if (is_vcf) {
		puts("##fileformat=VCFv4.1");
		if (fai) {
			printf("##reference=%s\n", fname);
			int i, n = faidx_fetch_nseq(fai);
			for (i=0; i<n; i++) {
				const char *seq = faidx_iseq(fai,i);
				int len = faidx_seq_len(fai, seq);
				printf("##contig=<ID=%s,length=%d>\n", seq, len);
			}
		}
		puts("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		if (show_2strand) {
			puts("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
			puts("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
		} else puts("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
		fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", stdout);
		for (i = 0; i < n; ++i) printf("\t%s", argv[optind+i]);
		putchar('\n');
	}
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		if (bed && !bed_overlap(bed, h->target_name[tid], pos, pos + 1)) continue; // not overlapping BED
		for (i = aux.tot_dp = 0; i < n; ++i) aux.tot_dp += n_plp[i];
		if (last_tid != tid) {
			if (is_fa && last_tid >= 0)
				write_fa(&aux, h->target_name[last_tid], 0, max_dev, l_ref);
			if (fai) { // switch of chromosomes
				free(ref);
				ref = fai_fetch(fai, h->target_name[tid], &l_ref);
			}
			last_tid = tid; last_pos = -1; aux.len = 0;
		}
		if (aux.tot_dp) {
			int k, r = 15, shift = 0, qual;
			allele_t *a;
			if (aux.tot_dp + 1 > aux.max_dp) { // expand array
				aux.max_dp = aux.tot_dp + 1;
				kroundup32(aux.max_dp);
				aux.a = (allele_t*)realloc(aux.a, aux.max_dp * sizeof(allele_t));
			}
			a = aux.a;
			// collect alleles
			r = (ref && pos - beg < l_ref)? seq_nt16_table[(int)ref[pos - beg]] : 15; // the reference allele
			for (i = aux.n_a = 0; i < n; ++i)
				for (j = 0; j < n_plp[i]; ++j) {
					a[aux.n_a] = pileup2allele(&plp[i][j], baseQ, (uint64_t)i<<32 | j, r, trim_len);
					if (!a[aux.n_a].is_skip) ++aux.n_a;
				}
			if (aux.n_a == 0) continue; // no reads are good enough; zero effective coverage
			// drop alleles in the Lianti mode
			if (lt_mode) lt_drop_reads(aux.n_a, aux.a);
			// count alleles
			ks_introsort(allele, aux.n_a, aux.a);
			count_alleles(&aux, n, qual_as_depth);
			// squeeze out weak alleles
			for (i = k = 0; i < aux.n_a; ++i)
				if (aux.support[a[i].k] >= min_support)
					a[k++] = a[i];
			if (k < aux.n_a) {
				if (k == 0) continue; // no alleles are good enough
				aux.n_a = k;
				count_alleles(&aux, n, qual_as_depth);
			}

			if (var_only && aux.n_alleles == 1 && a[0].hash>>63 == 0) continue; // var_only mode, but no ALT allele; skip
			if (is_fa) { // FASTA output
				int del_supp, c, is_ambi = 0, sum_dp = 0;
				if (pos - beg >= aux.max_len) { // expand arrays
					aux.max_len = pos - beg + 1;
					kroundup32(aux.max_len);
					aux.seq = (char*)realloc(aux.seq, aux.max_len);
					aux.depth = (int*)realloc(aux.depth, aux.max_len * sizeof(int));
				}
				for (i = last_pos + 1; i < pos; ++i) // fill gaps
					aux.seq[i - beg] = 'n', aux.depth[i - beg] = 0;
				for (j = del_supp = 0; j < n_plp[0]; ++j) // count reads supporting a deletion at this position
					if (plp[0][j].is_del)
						del_supp += qual_as_depth? bam_get_qual(plp[0][j].b)[plp[0][j].qpos] : 1;
				if (majority_fa || rand_fa) {
					int allele;
					if (majority_fa) {
						int max = 0, max_k = -1, n_max = 0;
						for (k = max = 0; k < aux.n_alleles; ++k)
							if (aux.support[k] > max) max = aux.support[k], max_k = k;
						assert(max_k >= 0);
						for (k = n_max = 0; k < aux.n_alleles; ++k)
							if (aux.support[k] == max) ++n_max;
						if (n_max > 1) {
							int r;
							r = (int)(n_max * drand48());
							if (r == n_max) r = n_max - 1;
							for (k = n_max = 0; k < aux.n_alleles; ++k)
								if (aux.support[k] == max && n_max++ == r)
									max_k = k;
						}
						allele = max_k;
						if (del_supp > max) is_ambi = 1;
					} else {
						double r;
						int tot;
						for (k = tot = 0; k < aux.n_alleles; ++k) tot += aux.support[k];
						r = tot * drand48();
						for (k = tot = 0; k < aux.n_alleles && tot + aux.support[k] < r; ++k)
							tot += aux.support[k];
						allele = k < aux.n_alleles? k : aux.n_alleles - 1;
					}
					for (i = 0; i < aux.n_a; ++i)
						if (a[i].k == allele) break;
					assert(i < aux.n_a);
					c = a[i].b;
					if (c != 1 && c != 2 && c != 4 && c != 8) c = 15, is_ambi = 1;
				} else {
					if (del_supp >= min_support) is_ambi = 1;
					if (aux.n_alleles > 2) is_ambi = 1;
					for (i = c = 0; i < aux.n_a; ++i) c |= a[i].b;
				}
				for (i = 0; i < aux.n_a; ++i) sum_dp += qual_as_depth? a[i].q : 1;
				c = (r == 1 || r == 2 || r == 4 || r == 8) && c == r? char_x : seq_nt16_str[c];
				if (is_ambi) c = tolower(c);
				aux.seq[pos - beg] = c;
				aux.depth[pos - beg] = sum_dp;
				aux.len = pos - beg + 1;
			} else { // print VCF or allele summary
				fputs(h->target_name[tid], stdout); printf("\t%d", pos+1);
				if (is_vcf) {
					fputs("\t.\t", stdout);
					for (i = 0; i <= aux.max_del; ++i) // print the reference allele up to the longest deletion
						putchar(ref && pos + i < l_ref + beg? ref[pos + i - beg] : 'N');
					putchar('\t');
				} else printf("\t%c\t", ref && pos < l_ref + beg? ref[pos - beg] : 'N'); // print a single reference base
				// print alleles
				if (!is_vcf || a[0].hash>>63) { // print if there is no reference allele
					print_allele(&plp[a[0].pos>>32][(uint32_t)a[0].pos], l_ref, ref, pos - beg, aux.max_del, is_vcf);
					if (aux.n_alleles > 1) putchar(',');
				}
				for (i = k = 1; i < aux.n_a; ++i)
					if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) {
						print_allele(&plp[a[i].pos>>32][(uint32_t)a[i].pos], l_ref, ref, pos - beg, aux.max_del, is_vcf);
						if (++k != aux.n_alleles) putchar(',');
					}
				if (is_vcf && aux.n_alleles == 1 && a[0].hash>>63 == 0) putchar('.'); // print placeholder if there is only the reference allele
				// compute and print qual
				for (i = !(a[0].hash>>63), qual = 0; i < aux.n_alleles; ++i)
					qual = qual > aux.support[i]? qual : aux.support[i];
				if (is_vcf) printf("\t%d\t.\t.\tGT:%s", qual, show_2strand? "ADF:ADR" : "AD");
				// print counts
				shift = (is_vcf && a[0].hash>>63); // in VCF, if there is no ref allele, we need to shift the allele number
				for (i = k = 0; i < n; ++i, k += aux.n_alleles) {
					int max1 = 0, max2 = 0, a1 = -1, a2 = -1, *sum_q = &aux.cnt_supp[k];
					// estimate genotype
					for (j = 0; j < aux.n_alleles; ++j)
						if (sum_q[j] > max1) max2 = max1, a2 = a1, max1 = sum_q[j], a1 = j;
						else if (sum_q[j] > max2) max2 = sum_q[j], a2 = j;
					if (max1 == 0 || (min_support > 0 && max1 < min_support)) a1 = a2 = -1;
					else if (max2 == 0 || (min_support > 0 && max2 < min_support)) a2 = a1;
					// print genotypes
					if (a1 < 0) printf("\t./.:");
					else printf("\t%d/%d:", a1 + shift, a2 + shift);
					// print counts
					if (show_2strand) {
						if (shift) fputs("0,", stdout);
						for (j = 0; j < aux.n_alleles; ++j) {
							if (j) putchar(',');
							printf("%d", aux.cnt_strand[(k+j)<<1]);
						}
						putchar(':');
						if (shift) fputs("0,", stdout);
						for (j = 0; j < aux.n_alleles; ++j) {
							if (j) putchar(',');
							printf("%d", aux.cnt_strand[(k+j)<<1|1]);
						}
					} else {
						if (shift) fputs("0,", stdout);
						for (j = 0; j < aux.n_alleles; ++j) {
							if (j) putchar(',');
							printf("%d", aux.cnt_supp[k+j]);
						}
					}
				} // ~for(i)
				putchar('\n');
			} // ~else if(is_fa)
			last_pos = pos;
		} // ~if(aux.tot_dp)
	} // ~while()
	if (is_fa && last_tid >= 0)
		write_fa(&aux, h->target_name[last_tid], 0, max_dev, l_ref);
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

	bam_hdr_destroy(h);
	for (i = 0; i < n; ++i) {
		bgzf_close(data[i]->fp);
		if (data[i]->itr) bam_itr_destroy(data[i]->itr);
		free(data[i]);
	}
	if (ref) free(ref);
	if (fai) fai_destroy(fai);
	free(aux.cnt_strand); free(aux.cnt_supp); free(aux.a);
	free(aux.seq); free(aux.depth);
	free(data); free(reg);
	if (bed) bed_destroy(bed);
	return 0;
}
