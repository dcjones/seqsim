
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
#include <algorithm>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>
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
    double       size_mean;
    double       size_std;
    double       nonunif_var;
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


/* A sparse vector giving the non-uniformity across the genome. */
class bias
{
    public:
        bias(const deque<transcript>& T)
        {
            printf("generating non-uniformity ... "); fflush(stdout);
            map<string, deque<exon> > all_exons;

            deque<transcript>::const_iterator t;
            for (t = T.begin(); t != T.end(); ++t) {
                all_exons[t->seqname].insert(all_exons[t->seqname].end(),
                        t->exons.begin(), t->exons.end());
            }

            map<string, deque<exon> >::iterator i;
            list<exon>::iterator u, v, w;
            for (i = all_exons.begin(); i != all_exons.end(); ++i) {
                sort(i->second.begin(), i->second.end());

                /* create a list, so we can efficiently erase redundant exons */
                list<exon> exons(i->second.begin(), i->second.end());

                u = exons.begin();
                v = exons.begin();
                ++v;

                while (v != exons.end()) {

                    while (v != exons.end() &&
                           u->start <= v->end && v->start <= u->end) {

                        u->end = max(u->end, v->end);
                        ++v;
                    }

                    // erase
                    w = v;
                    --w;

                    if (u != w) {
                        ++u;
                        exons.erase(u, v);
                    }

                    u = v;
                    ++v;
                }

                i->second.clear();
                i->second.insert(i->second.end(), exons.begin(), exons.end());
            }


            deque<exon>::const_iterator j;
            size_t k;
            for (i = all_exons.begin(); i != all_exons.end(); ++i) {
                for (j = i->second.begin(); j != i->second.end(); ++j) {
                    indexes[i->first].push_back(j->start);
                    blocks[i->first].push_back(vector<double>(j->end - j->start + 1));

                    /* entirely intependent bias: draw random beta variates */
                    vector<double>& bs = blocks[i->first].back();
                    if (args.nonunif_var > 0.0) {
                        double v = args.nonunif_var * (args.frag_pr * (1.0 - args.frag_pr));

                        double alpha = (args.frag_pr * (1.0 - args.frag_pr)) / v - 1.0;
                        double beta  = (1.0 - args.frag_pr) * alpha;
                        alpha = args.frag_pr * alpha;

                        assert(alpha > 0.0);
                        assert(beta  > 0.0);

                        for (k = 0; k < bs.size(); ++k) {
                            bs[k] = gsl_ran_beta(args.rng, alpha, beta);
                        }
                    }
                    else {
                        for (k = 0; k < bs.size(); ++k) {
                            bs[k] = args.frag_pr;
                        }
                    }
                }
            }

            printf("done.\n");
        }


        /* Build the cumulative sum of transcript t's bias vector. */
        void get_bias(vector<double>& bs, const transcript& t)
        {
            bs.resize(t.exonic_length());

            size_t i = 0;
            size_t j;
            size_t k;
            set<exon>::const_iterator e;
            for (e = t.exons.begin(); e != t.exons.end(); ++e) {

                /* 1. find the proper block in 'blocks' */
                j = find_block_idx(t.seqname, e->start);


                assert(indexes[t.seqname][j] <= e->start);

                k = e->start - indexes[t.seqname][j];

                block& B = blocks[t.seqname][j];

                assert((pos_t) (B.size() - k) >=
                        e->end - e->start);

                /* 2. copy the proper part of this block */
                copy(B.begin() + k,
                     B.begin() + k + (e->end - e->start + 1),
                     bs.begin() + i);

                i += e->end - e->start + 1;
            }
        }


    private:
        int find_block_idx(const string& seqname, pos_t start)
        {
            index_deque& I = indexes[seqname];
            index_deque::iterator i;

            i = lower_bound(I.begin(), I.end(), start);

            if (*i != start) --i;

            return (int) (i - I.begin());
        }

        typedef vector<double> block;
        typedef deque<block>   block_deque;
        typedef deque<int>     index_deque;

        map<string, block_deque> blocks;
        map<string, index_deque> indexes;
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


/* Generate a fragment count for each trancript, given m fragment in total. */
void generate_fragment_counts(unsigned int* cs,
                              const deque<transcript>& T,
                              const double* xs,
                              size_t m)
{
    printf("generating fragment counts ... "); fflush(stdout);

    const size_t n = T.size();

    /* probability of drawing a fragment from each transcript */
    double* ws = new double [n];

    size_t i;
    pos_t len, k;
    for (i = 0; i < n; ++i) {
        ws[i] = 0.0;
        len = T[i].exonic_length();

        const pos_t k_max = min(len, (pos_t) (args.size_mean + 6.0 * args.size_std));

        for (k = 1; k <= k_max; ++k) {
            ws[i] += (double) (len - k + 1) * gsl_ran_gaussian_pdf((double) k - args.size_mean, args.size_std);
        }

        ws[i] *= xs[i];
    }


    memset(cs, 0, n * sizeof(unsigned int));
    gsl_ran_multinomial(args.rng, n, (unsigned int) m, ws, cs);

    delete [] ws;

    printf("done.\n");
}



void generate_fragments(int* fs,
                        const deque<transcript>& T,
                        const unsigned int* cs)
{
    printf("generating fragments ...\n"); fflush(stdout);

    bias B(T);
    vector<double> as;
    vector<double> bs;

    const size_t n = T.size();
    size_t i, j, k;
    double a, p;
    double M = 1.0 / gsl_ran_gaussian_pdf(0.0, args.size_std);

    vector<double>::iterator l;

    pos_t start, end;

    k = 0;
    for (i = 0; i < n; ++i) {
        B.get_bias(as, T[i]);
        bs = as;

        // cumulative sum, for binary search to choose start position
        for (j = 1; j < as.size(); ++j) {
            as[j] += as[j - 1];
        }

        // cumulative sum of log(1-p), for binary searches for end position
        for (j = 0; j < bs.size(); ++j) {
            bs[j] = log(1.0 - bs[j]);
        }

        for (j = 1; j < bs.size(); ++j) {
            bs[j] += bs[j - 1];
        }

        for (j = 0; j < cs[i]; ) {
            /* choose a random start position */
            a = gsl_rng_uniform(args.rng) * as.back();
            l = lower_bound(as.begin(), as.end(), a);
            start = l - as.begin();

            /* choose a random end position */
            a = log(gsl_rng_uniform(args.rng)) + bs[start];
            l = lower_bound(bs.begin(), bs.end(), a, greater<double>());
            if (l == bs.end()) continue;
            end = l - bs.begin();

            /* simulate size selection: reject fragments according to a normal
             * distribution */
            p = M * gsl_ran_gaussian_pdf((double) (end - start + 1) - args.size_mean,
                                         args.size_std);
            if (gsl_rng_uniform(args.rng) > p) continue;

            fs[k * 3 + 0] = i;
            fs[k * 3 + 1] = start;
            fs[k * 3 + 2] = end;

            ++j;
            ++k;

            if (k % 100000 == 0) {
                printf("\t%zu\n", k);
            }
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

        if (start < l) {
            a = exon->start + start;
            b = exon->start + min(end + 1, l);

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
            "      --nonunif-std=N degree of nonuniformity in fragment generation\n"
            "\n");
}


int main(int argc, char* argv[])
{
    args.out         = "seqsim";
    args.n           = 25000000;
    args.paired_end  = false;
    args.stranded    = false;
    args.readlen     = 75;
    args.expr_pr     = 0.4;
    args.expr_a      = 1.0;
    args.frag_pr     = 0.004;
    args.frag_n      = 500000;
    args.size_mean   = 200.0;
    args.size_std    = 20.0;
    args.nonunif_var = 0.0;
    args.genome_fn   = NULL;
    args.genes_fn    = NULL;
    args.rng         = gsl_rng_alloc(gsl_rng_mt19937);

    setbuf(stdout, NULL);

    static struct option long_options[] =
    {
        {"help",        no_argument,       NULL, 'h'},
        {"out",         required_argument, NULL, 'o'},
        {"paired-end",  no_argument,       NULL, 'p'},
        {"stranded",    no_argument,       NULL, 's'},
        {"readlen",     required_argument, NULL, 'k'},
        {"expr-pr",     required_argument, NULL,  0 }, // 5
        {"expr-a",      required_argument, NULL,  0 },
        {"frag-pr",     required_argument, NULL,  0 },
        {"frag-n",      required_argument, NULL,  0 },
        {"size-mean",   required_argument, NULL,  0 },
        {"size-std",    required_argument, NULL,  0 },
        {"nonunif-var", required_argument, NULL,  0 },
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
                        args.size_mean = atof(optarg);
                        break;

                    case 10:
                        args.size_std = atof(optarg);
                        break;

                    case 11:
                        args.nonunif_var = atof(optarg);
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

    /* generate fragment counts */
    unsigned int* cs = new unsigned int [n];
    generate_fragment_counts(cs, T, xs, args.frag_n);

    /* generate random fragments */
    int* fs = new int [3 * args.frag_n];
    generate_fragments(fs, T, cs);

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






