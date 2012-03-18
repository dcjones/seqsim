
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
        "  -h, --help         print this message\n"
        "  -o, --output=FILE  output to the given file (output is to standard output by default)\n"
        "  -s, --seed=SEED    use the given seed to the random number generator.\n"
        );
}


static double frag_start_pr(const vector<double>& F, size_t len)
{
    if (len >= F.size()) len = F.size() - 1;

    double p = 0.0;
    size_t i;
    for (i = 0; i < len; ++i) p += F[i];
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


static void print_read(FILE* fout,
                       const transcript& t, const char* src,
                       pos_t start, pos_t end, uint8_t strand,
                       int mate)
{
    static uint64_t readnum = 0;

    if (mate != 2) ++readnum;

    fprintf(fout, ">seqsim.%"PRIu64"\n", readnum);

    /* TODO print sequence */



    fprintf(fout, "+\n");
    pos_t i;
    for (i = start; i <= end; ++i) fputc(max_qual, fout);
    fputc('\n', fout);
}


static void generate_fragment(
    gsl_rng* rng, params& P, FILE* fout,

    /* transcript from which the fragment originates */
    const transcript& t,

    /* cumulative distribution over fragment lengths */
    const vector<double>& F,

    /* cumulative distribution of fragment start positions */
    const vector<double>& fs,

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
    r = rmax / gsl_rng_uniform(rng);

    i_ = lower_bound(F.begin(), F.end(), r);
    size_t j = i + (i_ - F.begin());


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
            print_read(fout, t, seq, i, i + P.readlen - 1, strand_pos, 1);
            print_read(fout, t, seq, i, j - P.readlen + 1, strand_neg, 2);
        }
        else {
            print_read(fout, t, seq, i, j - P.readlen + 1, strand_neg, 1);
            print_read(fout, t, seq, i, i + P.readlen - 1, strand_pos, 2);
        }
    }
    else {
        if (strand == strand_pos) {
            print_read(fout, t, seq, i, i + P.readlen - 1, strand, 0);
        }
        else {
            print_read(fout, t, seq, j - P.readlen + 1, j, strand, 0);
        }
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

    for (i = 0; i < n; ++i) {
        a = gsl_cdf_gaussian_P(
                (double) (i - P.size_upper),
                P.size_sd);

        b = gsl_cdf_gaussian_P(
                (double) (i - P.size_lower),
                P.size_sd);

        F[i] = b - a;
    }

    return F;
}


int seqsim_generate(int argc, char* argv[])
{
    static struct option long_options[] =
    {
        {"help",   no_argument,       NULL, 'h'},
        {"output", required_argument, NULL, 'o'},
        {"seed",   required_argument, NULL, 's'},
        {0, 0, 0, 0}
    };

    unsigned long rng_seed = 135792468;
    char* out_fn = NULL;

    int opt, opt_idx;

    while (true) {
        opt = getopt_long(argc, argv, "ho:s:", long_options, &opt_idx);
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
    while (fgets(line, sizeof(line), expr_f)) {
        if ((sep1 = strchr(line, '\t')) == NULL) continue;
        if ((sep2 = strchr(sep1 + 1, '\t')) == NULL) continue;

        /* TODO: what form do we need this in? */
    }
    fclose(expr_f);


    /* Generate fragment length distribution. */
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
        z += (R[i] = effective_transcript_length(F, T[i].exonic_length()));
    }

    for (i = 0; i < n; ++i) R[i] /= z;


    /* Generate fragment counts for each transcript. */
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    vector<unsigned int> C;
    C.resize(n);
    gsl_ran_multinomial(rng, P.N, n, &R.at(0), &C.at(0));


    /* Index transcripts by chromosome. */
    map<string, vector<unsigned int> > trans_seqname;
    for (i = 0; i < n; ++i) {
        trans_seqname[T[i].seqname].push_back(i);
    }


    /* Generate fragments from each transcript. */
    fastq_t* fqf = fastq_open(genome_f);
    seq_t* seq = fastq_alloc_seq();

    /* transcript specific distribution over fragment start positions. */
    vector<double> fs;

    while (fastq_next(fqf, seq)) {
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

            while (C[*t]--) {
                generate_fragment(rng, P, stdout, T[*t], fraglen_cumdist, fs, seq->seq.s);
            }
        }
    }

    /* Print reads whose sequence is not in the genome for whatever reason. */
    // TODO

    fastq_close(fqf);
    fastq_free_seq(seq);

    fclose(genome_f);
    gsl_rng_free(rng);

    return EXIT_SUCCESS;
}

