
#include "generate.hpp"
#include "params.hpp"
#include "transcripts.hpp"
#include "fastq-parse.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cstdlib>
#include <cctype>
#include <inttypes.h>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;


/* All qualities are set to the highest on the phred=33 scale
 * as implemented in the Illumina pipeline. */
const char max_qual = 'G';

/* Min quality is used in places where there are 'N's in the sequence. */
const char min_qual = '!';


void seqsim_generate_usage(FILE* fout)
{
    fprintf(fout,
        "Usage: seqsim generate [options] <model.yaml> <genome.fa> <genes.gtf> <expression.tab>\n");
}


void seqsim_generate_help(FILE* fout)
{
    seqsim_generate_usage(fout);
    fprintf(fout,
        "Generate simulated RNA-Seq reads.\n\n"
        "Options:\n"
        "  -p, --prefix=PRE   output files should have this prefix (default: \"seqsim\")\n"
        "  -h, --help         print this message\n"
        "  -o, --output=FILE  output to the given file (output is to standard output by default)\n"
        "  -s, --seed=SEED    use the given seed to the random number generator.\n"
        );
}


static char complement(char c)
{
    switch( c ) {
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'n': return 'n';
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default:  return 'n';
    }
}

void seqrc(char* seq, int n)
{
    char c;
    int i,j;
    i = 0;
    j = n-1;
    while( i < j ) {
        c = complement(seq[i]);
        seq[i] = complement(seq[j]);
        seq[j] = c;
        i++; j--;
    }

    if( i == j ) seq[i] = complement(seq[i]);
}



static double frag_start_pr(const vector<double>& F, size_t len)
{
    static map<size_t, double> memo;

    if (len >= F.size()) return 0.0;

    map<size_t, double>::iterator ans = memo.find(len);
    if (ans != memo.end()) return ans->second;

    double p = 0.0;
    size_t i;
    for (i = 0; i < len; ++i) p += F[i];

    memo[len] = p;

    return p;
}


static void frag_start_cumdist(const vector<double>& F, vector<double>& fs)
{
    double z = 0.0;
    size_t len = fs.size();
    size_t i;
    for (i = 0; i < len; ++i) {
        z += (fs[i] = frag_start_pr(F, len - i));
    }

    for (i = 0; i < len; ++i) fs[i] /= z;
    for (i = 1; i < len; ++i) fs[i] += fs[i - 1];
}


static double effective_transcript_length(const vector<double>& F, size_t len)
{
    double p = 0.0;
    size_t i;
    for (i = 0; i < len; ++i) p += frag_start_pr(F, len - i);
    return p;
}


static void get_read_seq(char* dest,
                         const transcript& t, const char* src,
                         pos_t start, pos_t end, uint8_t strand)
{
    pos_t readlen = end - start + 1;

    pos_t u = 0; // offset into dest
    pos_t v = 0; // offset into src
    pos_t l;
    pos_t a, b;
    set<exon>::iterator e;
    for (e = t.exons.begin(); e != t.exons.end(); ++e) {
        l = e->end - e->start + 1;

        if (start < l) {
            a = e->start + start;
            b = e->start + min(end + 1, l);

            copy(src + a, src + b, dest + u);
            u += b - a;
            start += b - a;
        }

        /* make (start, end) relative to the next exon */
        start -= l;
        end   -= l;

        v += l;
        if (start < 0 || end < 0) break;
    }

    dest[u] = '\0';

    assert(strlen(dest) == (size_t) readlen);

    if (strand == strand_neg) seqrc(dest, readlen);
}


static void print_quals(FILE* fout, char* seq)
{
    while (*seq) {
        fputc(*seq == 'N' ? min_qual : max_qual, fout);
        ++seq;
    }
    fputc('\n', fout);
}



void generate_null_fragment(
    const params& P, FILE* fout1, FILE* fout2, uint64_t readnum)
{
    size_t i;
    if (P.paired) {
        fprintf(fout1, "@seqsim.%"PRIu64"/1\n", readnum);
        for (i = 0; i < P.readlen; ++i) fputc('N', fout1);
        fputs("\n+\n", fout1);
        for (i = 0; i < P.readlen; ++i) fputc('!', fout1);
        fputc('\n', fout1);


        fprintf(fout2, "@seqsim.%"PRIu64"/2\n", readnum);
        for (i = 0; i < P.readlen; ++i) fputc('N', fout2);
        fputs("\n+\n", fout2);
        for (i = 0; i < P.readlen; ++i) fputc('!', fout2);
        fputc('\n', fout2);
    }
    else {
        fprintf(fout1, "@seqsim.%"PRIu64"\n", readnum);
        for (i = 0; i < P.readlen; ++i) fputc('N', fout1);
        fputs("\n+\n", fout1);
        for (i = 0; i < P.readlen; ++i) fputc('!', fout1);
        fputc('\n', fout1);
    }
}



/* Generate a fragment from the given transcript. */
static void generate_fragment(
    gsl_rng* rng, params& P,

    /* output fastq files */
    FILE* fout1,
    FILE* fout2,

    uint64_t readnum,

    /* transcript from which the fragment originates */
    const transcript& t,

    /* cumulative distribution over fragment lengths */
    const vector<double>& F,

    /* cumulative distribution of fragment start positions */
    const vector<double>& fs,

    /* space for read sequences */
    char* readseq1,
    char* readseq2,

    /* nucleotide sequences from which the transcript originates */
    const char* seq)
{
    assert((size_t) t.exonic_length() == fs.size());
    size_t len = fs.size();


    /* choose start position (by binary search in fs) */
    double r = gsl_rng_uniform(rng);

    vector<double>::const_iterator i_ =
        lower_bound(fs.begin(), fs.end(), r);
    size_t i = i_ - fs.begin();


    /* choose end position */
    double rmax = F[len - i < F.size() ? len - i : F.size() - 1];
    r = rmax * gsl_rng_uniform(rng);

    i_ = lower_bound(F.begin(), F.end(), r);
    size_t j = i + (i_ - F.begin());

    /* nudge if the length is too far due to floating point imprecision */
    if (j >= len) j = len   - 1;


    /* choose strand */
    strand_t strand;
    if (gsl_ran_bernoulli(rng, P.strand_specificity)) {
        strand = t.strand;
    }
    else {
        strand = other_strand(t.strand);
    }


    /* extract sequence */
    if (P.paired) {
        if (strand == strand_pos) {
            get_read_seq(readseq1, t, seq, i, i + P.readlen - 1, strand_pos);
            get_read_seq(readseq2, t, seq, j - P.readlen + 1, j, strand_neg);
        }
        else {
            get_read_seq(readseq1, t, seq, j - P.readlen + 1, j, strand_neg);
            get_read_seq(readseq2, t, seq, i, i + P.readlen - 1, strand_pos);
        }

        fprintf(fout1, "@seqsim.%"PRIu64"/1\n%s\n+\n", readnum, readseq1);
        print_quals(fout1, readseq1);

        fprintf(fout2, "@seqsim.%"PRIu64"/2\n%s\n+\n", readnum, readseq2);
        print_quals(fout2, readseq2);
    }
    else {
        if (strand == strand_pos) {
            get_read_seq(readseq1, t, seq, i, i + P.readlen - 1, strand);
        }
        else {
            get_read_seq(readseq1, t, seq, j - P.readlen + 1, j, strand);
        }

        fprintf(fout1, "@seqsim.%"PRIu64"\n%s\n+\n", readnum, readseq1);
        print_quals(fout1, readseq1);
    }
}


/* Generate an explicit fragment length distribution. */
vector<double> make_frag_len_dist(const params& P)
{
    int n = P.size_upper + 10.0 * P.size_sd;
    vector<double> F;
    F.resize(n);
    int i;
    double a, b;
    double z = 0.0;

    for (i = 0; i < n; ++i) {
        a = gsl_cdf_gaussian_P(
                (double) (i - P.size_upper),
                P.size_sd);

        b = gsl_cdf_gaussian_P(
                (double) (i - P.size_lower),
                P.size_sd);

        z += (F[i] = b - a);
    }

    for (i = 0; i < n; ++i) F[i] /= z;

    return F;
}


int seqsim_generate(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",   no_argument,       NULL, 'h'},
        {"output", required_argument, NULL, 'o'},
        {"seed",   required_argument, NULL, 's'},
        {"prefix", required_argument, NULL, 'p'},
        {0, 0, 0, 0}
    };

    const char* out_prefix = "seqsim";
    unsigned long rng_seed = 135792468;
    char* out_fn = NULL;

    int opt, opt_idx;

    while (true) {
        opt = getopt_long(argc, argv, "ho:s:p:", long_options, &opt_idx);
        if (opt == -1) break;

        switch (opt) {
            case 'h':
                seqsim_generate_help(stdout);                
                return EXIT_SUCCESS;
                break;

            case 'o':
                out_fn = optarg;
                break;

            case 's':
                rng_seed = strtoul(optarg, NULL, 10);
                break;

            case 'p':
                out_prefix = optarg;
                break;

            default:
                abort();
        }
    }

    ++opt_idx;

    if (argc - opt_idx < 4) {
        seqsim_generate_usage(stderr);
        return EXIT_FAILURE;
    }

    params P;
    P.read(argv[opt_idx++]);

    /* Genome */
    FILE* genome_f = fopen(argv[opt_idx], "r");
    if (genome_f == NULL) {
        fprintf(stderr, "seqsim generate: can't open FASTA file \"%s\".", argv[opt_idx]);
        exit(EXIT_FAILURE);
    }
    opt_idx++;


    /* Genes */
    FILE* genes_f = fopen(argv[opt_idx], "r");
    if (genes_f == NULL) {
        fprintf(stderr, "seqsim generate: can't open GTF file \"%s\".", argv[opt_idx]);
        exit(EXIT_FAILURE);
    }

    vector<transcript> T;
    read_transcripts_from_gtf(genes_f, T);
    fclose(genes_f);
    size_t n = T.size();
    opt_idx++;


    /* Expression. */
    FILE* expr_f = fopen(argv[opt_idx], "r");
    if (expr_f == NULL) {
        fprintf(stderr, "seqsim generate: can't open expression file \"%s\".", argv[opt_idx]);
    }

    char  line[512];
    char *sep1, *sep2;
    map<string, double> expr;
    while (fgets(line, sizeof(line), expr_f)) {
        if ((sep1 = strchr(line, '\t')) == NULL) continue;
        if ((sep2 = strchr(sep1 + 1, '\t')) == NULL) continue;

        string tid = string(sep1 + 1, sep2 - sep1 - 1);
        double e = atof(sep2 + 1);
        expr[tid] = e;
    }
    fclose(expr_f);


    /* Generate fragment length distribution. */
    fprintf(stderr, "calculating sampling rate ... ");
    vector<double> F = make_frag_len_dist(P);
    vector<double> fraglen_cumdist;
    fraglen_cumdist.resize(F.size());
    size_t i;
    fraglen_cumdist[0] = F[0];
    for (i = 1; i < F.size(); ++i) {
        fraglen_cumdist[i] = F[i] + fraglen_cumdist[i - 1];
    }


    /* Compute a sampling rate for each transcript, proportional
     * to the probability of drawing a read from that transcript.
     */
    vector<double> R;
    R.resize(n);

    double z = 0.0;
    for (i = 0; i < n; ++i) {
        R[i] = effective_transcript_length(F, T[i].exonic_length());
        R[i] *= expr[T[i].transcript_id];
        z += R[i];
    }

    for (i = 0; i < n; ++i) R[i] /= z;


    /* Generate fragment counts for each transcript. */
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    vector<unsigned int> C;
    C.resize(n);
    gsl_ran_multinomial(rng, n, P.N, &R.at(0), &C.at(0));


    /* output fragment counts */
    string fn = out_prefix;
    fn += ".frag_count.tab";
    FILE* frag_cnt_f = fopen(fn.c_str(), "w");
    if (frag_cnt_f == NULL) {
        fprintf(stderr,
            "seqsim generate: can't open file \"%s\" for writing.\n", fn.c_str());
        return EXIT_FAILURE;
    }

    for (i = 0; i < n; ++i) {
        fprintf(frag_cnt_f, "%s\t%u\n",
            T[i].transcript_id.c_str(), C[i]);
    }

    fclose(frag_cnt_f);

    fprintf(stderr, "done.\n");


    /* Index transcripts by chromosome. */
    map<string, vector<unsigned int> > trans_seqname;
    for (i = 0; i < n; ++i) {
        trans_seqname[T[i].seqname].push_back(i);
    }


    /* output files */
    FILE *fout1, *fout2;

    if (P.paired) {
        string fn = out_prefix;
        fn += "_1.fastq";
        fout1 = fopen(fn.c_str(), "w");
        if (fout1 == NULL) {
            fprintf(stderr,
                "seqsim generate: can't open file \"%s\" for writing.\n", fn.c_str());
            return EXIT_FAILURE;
        }

        fn = out_prefix;
        fn += "_2.fastq";
        fout2 = fopen(fn.c_str(), "w");
        if (fout2 == NULL) {
            fprintf(stderr,
                "seqsim generate: can't open file \"%s\" for writing.\n", fn.c_str());
            return EXIT_FAILURE;
        }
    }
    else {
        string fn = out_prefix;
        fn += ".fastq";
        fout1 = fopen(fn.c_str(), "w");
        if (fout1 == NULL) {
            fprintf(stderr,
                "seqsim generate: can't open file \"%s\" for writing.\n", fn.c_str());
            return EXIT_FAILURE;
        }

        fout2 = NULL;
    }

    /* space for read sequences */
    char* readseq1 = new char[P.readlen];
    char* readseq2 = new char[P.readlen];


    /* Generate fragments from each transcript. */
    fastq_t* fqf = fastq_open(genome_f);
    seq_t* seq = fastq_alloc_seq();

    /* transcript specific distribution over fragment start positions. */
    vector<double> fs;
    uint64_t readnum = 0;

    fprintf(stderr, "generating reads ...\n");

    while (fastq_next(fqf, seq)) {
        fprintf(stderr, "\t%s ...\n", seq->id1.s);

        map<string, vector<unsigned int> >::iterator j = trans_seqname.find(seq->id1.s);
        if (j == trans_seqname.end()) continue;
        vector<unsigned int>& ts = j->second;


        char* c = seq->seq.s;
        while (*c) {
            *c = toupper(*c);
            c++;
        }


        vector<unsigned int>::iterator t;
        for (t = ts.begin(); t != ts.end(); ++t) {
            if (C[*t] > 0) {
                fs.resize(T[*t].exonic_length());
                frag_start_cumdist(F, fs);
            }

            while (C[*t] > 0) {
                --C[*t];

                generate_fragment(
                    rng, P, fout1, fout2,
                    ++readnum,
                    T[*t], fraglen_cumdist, fs,
                    readseq1, readseq2,
                    seq->seq.s);
            }
        }
    }

    /* Print reads whose sequence is not in the genome for whatever reason. */
    for (i = 0; i < n; ++i) {
        while (C[i] > 0) {
            --C[i];
            generate_null_fragment(P, fout1, fout2, ++readnum);
        }
    }

    if (fout1) fclose(fout1);
    if (fout2) fclose(fout2);

    fastq_close(fqf);
    fastq_free_seq(seq);

    fclose(genome_f);
    gsl_rng_free(rng);

    fprintf(stderr, "done.\n");

    return EXIT_SUCCESS;
}

