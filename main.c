#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

#define LT_VERSION "r84"

int main_trim(int argc, char *argv[]);
int main_ldup(int argc, char *argv[]);
int main_group(int argc, char *argv[]);
int main_count(int argc, char *argv[]);
int main_cnv(int argc, char *argv[]);
int main_pileup(int argc, char *argv[]);
int main_lorenz(int argc, char *argv[]);

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int main(int argc, char *argv[])
{
	int ret = 0, i;
	double t_start;
	liftrlimit();
	if (argc == 1) {
		fprintf(stderr, "Usage: lianti <command> <arguments>\n");
		fprintf(stderr, "Commands:\n");
		fprintf(stderr, "  trim     trim binding motifs and adapter sequences\n");
		fprintf(stderr, "  ldup     mark Illumina PCR duplicates\n");
		fprintf(stderr, "  group    group reads into alleles\n");
		fprintf(stderr, "  count    compute allele depth\n");
		fprintf(stderr, "  cnv      call copy number variations\n");
		fprintf(stderr, "  pileup   lianti-aware pileup\n");
		fprintf(stderr, "  lorenz   compute the Lorenz evaluation curve\n");
		fprintf(stderr, "  version  print version number\n\n");
		fprintf(stderr, "Typical workflow:\n");
		fprintf(stderr, "  seqtk mergepe read1.fq.gz read2.fq.gz | lianti trim - | bwa mem -Cpt8 ref.fa - \\\n");
		fprintf(stderr, "    | samtools view -uS - | sambamba sort /dev/stdin | lianti ldup - > aln.bam\n");
		fprintf(stderr, "  lianti group aln.bam | bgzip > alleles.bed.gz\n");
		fprintf(stderr, "  lianti count alleles.bed.gz > depth.bed.gz\n");
		return 1;
	}
	t_start = realtime();
	if (strcmp(argv[1], "trim") == 0) ret = main_trim(argc-1, argv+1);
	else if (strcmp(argv[1], "ldup") == 0) ret = main_ldup(argc-1, argv+1);
	else if (strcmp(argv[1], "group") == 0) ret = main_group(argc-1, argv+1);
	else if (strcmp(argv[1], "count") == 0) ret = main_count(argc-1, argv+1);
	else if (strcmp(argv[1], "cnv") == 0) ret = main_cnv(argc-1, argv+1);
	else if (strcmp(argv[1], "pileup") == 0) ret = main_pileup(argc-1, argv+1);
	else if (strcmp(argv[1], "lorenz") == 0) ret = main_lorenz(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(LT_VERSION);
		return 0;
	} else {
		fprintf(stderr, "[E::%s] unknown command\n", __func__);
		return 1;
	}
	if (ret == 0) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, LT_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_start, cputime());
	}
	return ret;
}
