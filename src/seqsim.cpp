
/*
 * seqsim
 * ------
 *
 * This is an extremely simple RNA-Seq simulation. I make no guarantees of accuracy
 * or correctness. This is not a substitute for testing with real data!
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "config.h"
#include "common.h"
#include "gtf_parse.h"
#include "fastq-parse.h"
#include <deque>
#include <map>
#include <set>
#include <string>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <getopt.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

const char max_qual = 'g';


/* global arguments */
static struct
{
    const char*  out;
    size_t       n;
    bool         paired_end;
    bool         stranded;
    size_t       readlen;
    double       expr_pr;
    double       expr_a;
    double       frag_pr;
    size_t       frag_n;
    unsigned int size_low;
    unsigned int size_high;
    double       size_noise;
    const char*  genome_fn;
    const char*  genes_fn;
    gsl_rng*     rng;
} args;



/* Get all the sequence names from a FILE file. */
set<string> get_fasta_seqnames(const char* fn)
{
    printf("examining reference sequence ... ");

    set<string> seqnames;

    const size_t bufsiz = 4096;
    char* buf = new char [bufsiz];
    size_t n;
    
    FILE* fin = fopen_or_die(fn, "r");

    while (fgets(buf, bufsiz, fin)) {
        n = strlen(buf);
        if (buf[0] == '>') {
            buf[n - 1] = '\0';
            seqnames.insert(string(buf + 1));
        }
    }

    fclose(fin);
    delete [] buf;

    printf("done. (%zu sequences found)\n", seqnames.size());

    return seqnames;
}


/* Determine whether [u1, v1] overlaps [u2, v2] */
bool overlaps(pos_t u1, pos_t v1, pos_t u2, pos_t v2)
{
    assert(u1 <= v1);
    assert(u2 <= v2);

    return u1 <= v2 && v1 >= v2;
}


/* An interval, representing an exon. */
struct exon
{
    exon(pos_t start, pos_t end)
        : start(start), end(end)
    {
    }

    pos_t start;
    pos_t end;

    bool operator < (const exon& other) const
    {
        if (start != other.start) return start < other.start;
        else                      return end   < other.end;
    }
};


/* Representation of a single transcript. */
class transcript
{
    public:
        transcript()
            : strand(strand_na)
        {
        }

        /* Length of the mature mRNA. */
        pos_t exonic_length() const
        {
            pos_t len = 0;
            set<exon>::iterator i;
            for (i = exons.begin(); i != exons.end(); ++i) {
                len += i->end - i->start + 1;
            }

            return len;
        }

        /* Insert an exon into the transcript. */
        void add_exon(const gtf_row_t* row)
        {
            if (strcmp(row->feature->s, "exon") != 0) return;

            str_t* val;

            if (seqname.empty())     seqname = row->seqname->s;
            if (strand == strand_na) strand  = row->strand;
            if (gene_id.empty()) {
                val = (str_t*) str_map_get(row->attributes, "gene_id", 7);
                gene_id = val->s;
            }

            if (transcript_id.empty()) {
                val = (str_t*) str_map_get(row->attributes, "transcript_id", 13);
                transcript_id = val->s;
            }

            exons.insert(exon(row->start - 1, row->end - 1));
        }


        strand_t strand;
        string   seqname;
        string   gene_id;
        string   transcript_id;

        set<exon> exons;
};



/* Parse all the transcript from a GTF file. */
void read_transcripts(deque<transcript>& T, const char* fn, const set<string>& seqnames)
{
    printf("parsing gtf ... ");

    FILE* fin = fopen(fn, "r");

    gtf_file_t* gf = gtf_file_alloc(fin);
    gtf_row_t* row = gtf_row_alloc();

    set<string>::iterator j;

    map<string, transcript> ts;
    map<string, transcript>::iterator i;
    str_t* val;

    while (gtf_next(gf, row)) {
        if (seqnames.find(row->seqname->s) == seqnames.end()) continue;

        val = (str_t*) str_map_get(row->attributes, "transcript_id", 13);
        if (val && val->s) {
            ts[val->s].add_exon(row);
        }
    }

    gtf_row_free(row);
    gtf_file_free(gf);

    fclose(fin);

    for (i = ts.begin(); i != ts.end(); ++i) {
        T.push_back(i->second);
    }

    printf("done. (%zu transcripts)\n", T.size());
}


void generate_expression(double* xs, size_t n)
{
    printf("generating expression ... ");

    // choose the number of expressed genes
    size_t m = (size_t) gsl_ran_binomial(args.rng, args.expr_pr, n);

    // choose which transcripts are expressed
    int* idx = new int [n];
    size_t i;
    for (i = 0; i < n; ++i) idx[i] = i;
    gsl_ran_shuffle(args.rng, (void*) idx, n, sizeof(int));

    // generate expression values
    double* alpha = new double [m];
    for (i = 0; i < m; ++i) alpha[i] = args.expr_a;

    double* ys = new double [m];
    gsl_ran_dirichlet(args.rng, m, alpha, ys);

    memset(xs, 0, n * sizeof(double));
    for (i = 0; i < m; ++i) {
        xs[idx[i]] = ys[i];
    }

    delete [] ys;
    delete [] alpha;
    delete [] idx;

    printf("done.\n");
}


/* Print an expression table. */
void print_expression(const deque<transcript> T, const double* xs)
{
    char fn[1024];
    snprintf(fn, sizeof(fn), "%s.expr", args.out);

    FILE* fout = fopen_or_die(fn, "w");

    fprintf(fout, "gene_id\ttranscript_id\texpr\n");
    deque<transcript>::const_iterator t;
    size_t i;
    for (i = 0, t = T.begin(); t != T.end(); ++i, ++t) {
        fprintf(fout, "%s\t%s\t%0.6e\n",
                      t->gene_id.c_str(),
                      t->transcript_id.c_str(),
                      xs[i]);
    }

    fclose(fout);
}


/* binary search: return an index i, where for all j < i,
 *      xs[j] <= y
 *  and for all k > i
 *      xs[k] > y
 */
int search_sorted(double y, double* xs, size_t n)
{
    int a = 0;
    int b = (int) n - 1;
    int i = -1;

    while (a <= b) {
        i = a + (b - a) / 2;

        if (xs[i] < y)                   a = i + 1;
        else if(i > 0 && y <= xs[i - 1]) b = i - 1;
        else break;
    }

    return i;
}


void generate_fragments(int* fs,
                        const deque<transcript>& T,
                        const double* xs,
                        size_t m)
{
    printf("generating fragments ...\n");
    size_t n = T.size();

    // compute a weight vector
    double* ws = new double [n];
    size_t i;
    for (i = 0; i < n; ++i) {
        ws[i] = (double) T[i].exonic_length() * xs[i];
    }

    double z = 0;
    for (i = 0; i < n; ++i) z += ws[i];
    for (i = 0; i < n; ++i) ws[i] /= z;

    double cumsum = 0.0;
    z = 0.0;
    for (i = 0; i < n; ++i) {
        z = ws[i];
        ws[i] += cumsum;
        cumsum += z;
    }

    pos_t start, end;
    pos_t fraglen;
    pos_t l;

    double nudge;

    double r;
    size_t j;
    for (i = 0; i < m; ) {
        // choose a random transcript
        r = gsl_rng_uniform(args.rng);
        j = (size_t) search_sorted(r, ws, n);

        l = T[j].exonic_length();

        // choose start and end points
        start = gsl_rng_uniform_int(args.rng, l);
        end   = start + gsl_ran_geometric(args.rng, args.frag_pr);

        // reject improper fragments (end breakpoint outside the mRNA)
        if (end >= l) continue;

        // reject fragments too short to sequence
        fraglen = end - start + 1;
        if (fraglen < (pos_t) args.readlen) continue;

        // reject fragments lost in size selection
        nudge = gsl_ran_gaussian(args.rng, args.size_noise);
        if ((double) fraglen < (double) args.size_low + nudge) continue;

        nudge = gsl_ran_gaussian(args.rng, args.size_noise);
        if ((double) fraglen > (double) args.size_high + nudge) continue;

        fs[i * 3 + 0] = j;
        fs[i * 3 + 1] = start;
        fs[i * 3 + 2] = end;

        ++i;

        if (i % 100000 == 0) {
            printf("\t%zu\n", i);
        }
    }

    printf("done.\n");
}



/* Given n size-selected fragments, generate m reads. */
void generate_reads(int* rs, size_t n, size_t m)
{
    printf("generating reads ... ");

    memset(rs, 0, n * sizeof(int));

    size_t i;
    for (i = 0; i < m; ++i) {
        rs[gsl_rng_uniform_int(args.rng, n)] += 1;
    }

    printf("done.\n");
}



static void upper_str(char* s)
{
    while (*s != '\0') {
        *s = toupper(*s);
        ++s;
    }
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




void get_seq(char* dest, const char* src, pos_t start, pos_t end,
             const transcript& t)
{
    pos_t u = 0; // offset into dest
    pos_t v = 0; // offset into src
    pos_t l;

    pos_t a, b;

    set<exon>::iterator exon;
    for (exon = t.exons.begin(); exon != t.exons.end(); ++exon) {
        l = exon->end - exon->start + 1;

        if (v <= start && start < v + l) {
            a = exon->start + max(start, v);
            b = exon->start + min(end + 1, v + l);

            copy(src + a, src + b, dest + u);
            u += b - a;
            start = v + l;

            if (start > end) break;
        }

        v += l;
    }

    dest[u] = '\0';

    assert(strlen(dest) == args.readlen);
}



void print_read(FILE* fout1, FILE* fout2,
                size_t readnum, pos_t start, pos_t end,
                const transcript& t, const char* seq,
                char* seq1, char* seq2)
{
    get_seq(seq1, seq, start, start + args.readlen - 1, t);
    get_seq(seq2, seq, end - args.readlen + 1, end, t);

    upper_str(seq1);
    upper_str(seq2);

    strand_t strand;

    if (args.stranded) strand = t.strand;
    else {
        strand = gsl_rng_uniform(args.rng) < 0.5 ? strand_pos : strand_neg;
    }

    if (strand == strand_pos) {
        seqrc(seq2, args.readlen);
    }
    else {
        seqrc(seq2, args.readlen);
        swap(seq1, seq2);
    }


    if (args.paired_end) {
        fprintf(fout1, "@seqsim.%s.%zu/1\n%s\n", args.out, readnum, seq1);
        memset(seq1, max_qual, args.readlen);
        fprintf(fout1, "+\n%s\n", seq1);

        fprintf(fout2, "@seqsim.%s.%zu/2\n%s\n", args.out, readnum, seq2);
        memset(seq2, max_qual, args.readlen);
        fprintf(fout2, "+\n%s\n", seq2);
    }
    else {
        fprintf(fout1, "@seqsim.%s.%zu\n%s\n", args.out, readnum, seq1);
        memset(seq1, max_qual, args.readlen);
        fprintf(fout1, "+\n%s\n", seq1);
    }
}




void print_usage(FILE* fout)
{
    fprintf(fout, "Usage: seqsim [options] genome.fa genes.gtf\n");
}


void print_help(FILE* fout)
{
    print_usage(fout);
    fprintf(fout, "Version: %s\n\n", VERSION);
    fprintf(fout,
            "Options:\n"
            "  -h, --help          show this help message and exit\n"
            "  -n N                number of reads / pairs to generate\n"
            "  -o, --out=PRE       output files will have the given prefix\n"
            "                      (default: 'seqsim')\n"
            "  -s, --stranded      generate strand specific reads\n"
            "  -p, --paired-end    generate paired-end reads\n"
            "      --expr-pr=N     probability that a transcript is expressed\n"
            "      --expr_a=N      expression distribution parameter (hint: a\n"
            "                      small number results in more variance)\n"
            "      --frag-pr=N     probability of fragmentation\n"
            "      --frag-n=N      number of fragments to generate\n"
            "      --size-low=N    size selection lower bound\n"
            "      --size-high=N   size selection upper bound\n"
            "      --size-noise=N  inexectness / variance in size selection\n"
            "\n");
}


int main(int argc, char* argv[])
{
    args.out        = "seqsim";
    args.n          = 25000000;
    args.paired_end = false;
    args.stranded   = false;
    args.readlen    = 75;
    args.expr_pr    = 0.4;
    args.expr_a     = 1.0;
    args.frag_pr    = 0.004;
    //args.frag_n     = 50000000;
    args.frag_n     = 500000;
    args.size_low   = 180;
    args.size_high  = 220;
    args.size_noise = 5.0;
    args.genome_fn  = NULL;
    args.genes_fn   = NULL;
    args.rng        = gsl_rng_alloc(gsl_rng_mt19937);

    setbuf(stdout, NULL);

    static struct option long_options[] =
    {
        {"help",       no_argument,       NULL, 'h'},
        {"out",        required_argument, NULL, 'o'},
        {"paired-end", no_argument,       NULL, 'p'},
        {"stranded",   no_argument,       NULL, 's'},
        {"readlen",    required_argument, NULL, 'k'},
        {"expr-pr",    required_argument, NULL,  0 }, // 5
        {"expr-a",     required_argument, NULL,  0 },
        {"frag-pr",    required_argument, NULL,  0 },
        {"frag-n",     required_argument, NULL,  0 },
        {"size-low",   required_argument, NULL,  0 },
        {"size-high",  required_argument, NULL,  0 },
        {"size-noise", required_argument, NULL,  0 },
        {0, 0, 0, 0}
    };

    const char* short_options = "hn:opsk:";

    int opt;
    int opt_idx;

    while (1) {
        opt = getopt_long(argc, argv, short_options, long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 0:
                switch (opt_idx) {
                    case 5:
                        args.expr_pr = atof(optarg);
                        break;

                    case 6:
                        args.expr_a = atof(optarg);
                        break;

                    case 7:
                        args.frag_pr = atof(optarg);
                        break;

                    case 8:
                        args.frag_n = (size_t) atoi(optarg);
                        break;

                    case 9:
                        args.size_low = (unsigned int) atoi(optarg);
                        break;

                    case 10:
                        args.size_high = (unsigned int) atoi(optarg);
                        break;

                    case 11:
                        args.size_noise = atof(optarg);
                        break;
                };
                break;

            case 'h':
                print_help(stdout);
                return EXIT_SUCCESS;

            case 'n':
                args.n = atoi(optarg);
                break;

            case 'o':
                args.out = optarg;
                break;

            case 'p':
                args.paired_end = true;
                break;

            case 's':
                args.stranded = true;
                break;

            case 'k':
                args.readlen = atoi(optarg);
                break;

            case '?':
                return 1;

            default:
                abort();
        }
    }


    /* no positional arguments */
    if (optind == argc) {
        print_usage(stdout);
        return 0;
    }

    /* too few positional arguments */
    else if (optind + 1 == argc) {
        fprintf(stderr, "Too few arguments.\n\n");
        print_usage(stderr);
        return 1;
    }

    /* too many */
    else if (optind + 2 < argc) {
        fprintf(stderr, "Too many arguments.\n\n");
        print_usage(stderr);
        return 1;
    }

    args.genome_fn = argv[optind];
    args.genes_fn  = argv[optind + 1];

    set<string> seqnames = get_fasta_seqnames(args.genome_fn);

    deque<transcript> T;
    read_transcripts(T, args.genes_fn, seqnames);
    size_t n = T.size();

    /* generate / output true expression */
    double* xs = new double [n];
    generate_expression(xs, n);

    /* print expression */
    print_expression(T, xs);

    /* generate random fragments */
    int* fs = new int [3 * args.frag_n];
    generate_fragments(fs, T, xs, args.frag_n);

    /* sample reads from random fragments */
    int* rs = new int [args.frag_n];
    generate_reads(rs, args.frag_n, args.n);

    /* extract sequence */
    printf("getting sequence ...\n");

    FILE* fout1;
    FILE* fout2;
    char out_fn[1024];

    if (args.paired_end) {
        snprintf(out_fn, sizeof(out_fn), "%s.1.fastq", args.out);
        fout1 = fopen_or_die(out_fn, "w");

        snprintf(out_fn, sizeof(out_fn), "%s.2.fastq", args.out);
        fout2 = fopen_or_die(out_fn, "w");
    }
    else {
        snprintf(out_fn, sizeof(out_fn), "%s.fastq", args.out);
        fout1 = fopen_or_die(out_fn, "w");

        fout2 = NULL;
    }

    FILE* fin = fopen_or_die(args.genome_fn, "r");
    fastq_t* fqf = fastq_open(fin);


    seq_t* seq = fastq_alloc_seq();


    deque<transcript>::iterator t;
    int i, u, v;
    set<int> us;

    char* seq1 = new char [args.readlen + 1];
    char* seq2 = new char [args.readlen + 1];


    unsigned int readnum = 1;
    while (fastq_next(fqf, seq)) {
        printf("\t%s\n", seq->id1.s);

        // indices of transcripts from this sequence
        us.clear();
        for (u = 0, t = T.begin(); t != T.end(); ++u, ++t) {
            if (t->seqname == seq->id1.s) us.insert(u);
        }

        for (v = 0; v < (int) args.frag_n; ++v) {
            /* skip fragments from transcripts not on this sequence */
            if (us.find(fs[v * 3 + 0]) == us.end()) continue;

            /* print generate reads */
            for (i = 0; i < rs[v]; ++i, ++readnum) {
                print_read(fout1, fout2, readnum,
                           fs[v * 3 + 1], fs[v * 3 + 2],
                           T[fs[v * 3]], seq->seq.s,
                           seq1, seq2);
            }
        }
    }

    printf("done.\n");

    fclose(fout1);
    if (fout2) fclose(fout2);

    fastq_free_seq(seq);
    fastq_close(fqf);

    delete [] seq1;
    delete [] seq2;
    delete [] rs;
    delete [] fs;
    delete [] xs;
    gsl_rng_free(args.rng);

    return 0;
}






