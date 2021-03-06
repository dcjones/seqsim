
#include "express.hpp"
#include "params.hpp"
#include "transcripts.hpp"
#include "log_norm_mix.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <getopt.h>
#include <map>

using namespace std;

void seqsim_express_usage(FILE* fout)
{
    fprintf(fout, "Usage: seqsim express [options] <model.yaml> <genes.gtf>\n");
}


void seqsim_express_help(FILE* fout)
{
    seqsim_express_usage(fout);
    fprintf(fout,
        "Generate expression values for as set of genes and transcripts.\n\n"
        "Options:\n"
        "  -h, --help         print this message\n"
        "  -o, --output=FILE  output to the given file (output is to standard output by default)\n"
        "  -s, --seed=SEED    use the given seed to the random number generator.\n"
        );
}


int seqsim_express(int argc, char* argv[])
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

    int opt;

    while (true) {
        opt = getopt_long(argc, argv, "ho:s:", long_options, NULL);
        if (opt == -1) break;

        switch (opt) {
            case 'h':
                seqsim_express_help(stdout);
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

    // ++optind;

    if (argc - optind < 2) {
        seqsim_express_usage(stderr);
        return EXIT_SUCCESS;
    }

    params P;
    P.read(argv[optind]);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, rng_seed);

    FILE* gtf_fin = fopen(argv[optind + 1], "r");
    if (gtf_fin == NULL) {
        fprintf(stderr, "seqsim: %s: cannot open file.\n", argv[optind + 1]);
        return EXIT_FAILURE;
    }

    vector<transcript> T;
    read_transcripts_from_gtf(gtf_fin, T);
    fclose(gtf_fin);

    /* map gene IDs to transcript ids and transcript IDs to transcripts */
    map<string, set<unsigned int> > gids;
    vector<transcript>::iterator t;
    size_t i;
    for (t = T.begin(), i = 0; t != T.end(); ++t, ++i) {
        gids[t->gene_id].insert(i);
    }

    size_t n = T.size();    // number of transcripts
    size_t m = gids.size(); // number of genes

    fprintf(stderr, "generating expression values for %zu genes and %zu transcripts ... ",
                    m, n);


    /* generate gene-wise expression. */
    double z = 0.0;
    double* theta = new double [m];
    for (i = 0; i < m; ++i) {
        theta[i] = log_norm_mix(rng, P.gene_exp_k,
            &P.gene_exp_p.at(0), &P.gene_exp_mu.at(0), &P.gene_exp_sd.at(0));
        z += theta[i];
    }

    for (i = 0; i < m; ++i) theta[i] /= z;


    /* generate transcript-wise expression. */
    double* phi = new double [n];
    for (i = 0; i < n; ++i) {
        /* Note: drawing vector elements from gamma(a, 1) and normalizing
         * results in a symmetric dirichlet(a) distribution. */
        phi[i] = gsl_ran_gamma(rng, P.trans_exp_alpha, 1.0);
    }

    map<string, set<unsigned int> >::iterator gid;
    set<unsigned int>::iterator tid;

    for (i = 0, gid = gids.begin(); gid != gids.end(); ++i, ++gid) {
        z = 0.0;
        for (tid = gid->second.begin(); tid != gid->second.end(); ++tid) {
            z += phi[*tid];
        }

        for (tid = gid->second.begin(); tid != gid->second.end(); ++tid) {
            phi[*tid] = theta[i] * (phi[*tid] / z);
        }
    }

    /* print expression values */
    for (i = 0; i < n; ++i) {
        printf("%s\t%s\t%e\n",
            T[i].gene_id.c_str(),
            T[i].transcript_id.c_str(),
            phi[i]);
    }

    gsl_rng_free(rng);

    fprintf(stderr, "done.\n");

    return EXIT_SUCCESS;
}

