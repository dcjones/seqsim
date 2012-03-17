
#include "perturb.hpp"
#include "params.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <getopt.h>
#include <string>
#include <map>

using namespace std;


void seqsim_perturb_usage(FILE* fout)
{
    fprintf(fout, "Usage: seqsim perturb [options] <model.yml> <perturbation> <expression.tab>");
}

void seqsim_perturb_help(FILE* fout)
{
    seqsim_perturb_usage(fout);
    fprintf(fout,
        "Modify expression values according to particular model.\n\n"
        "Options:\n"
        "  -h, --help \n"
        "  -o, --output=FILE\n"
        "  -s, --seed=SEED\n"
        );
}



static void perturb_phi(
    gsl_rng* rng,
    const perturb_params& perturbation,
    const vector<string>& gene_ids,
    vector<double>& phi)
{
    size_t i;
    size_t n = gene_ids.size();

    vector<unsigned int> gene_idx;
    gene_idx.reserve(n);

    unsigned int idx;
    map<string, unsigned int> gene_id_set;
    for (i = 0; i < n; ++i) {
        if (gene_id_set.find(gene_ids[i]) == gene_id_set.end()) {
            idx = gene_id_set.size();
            gene_id_set[gene_ids[i]] = idx;
        }

        gene_idx.push_back(gene_id_set[gene_ids[i]]);
    }

    size_t m = gene_id_set.size();


    /* perturb gene expression */

    /* factor by which each gene's expression is modified */
    double* ds = new double[m];
    for (i = 0; i < m; ++i) {
        if (gsl_rng_uniform(rng) < perturbation.gene_pr) {
            ds[i] = gsl_ran_lognormal(
                rng, perturbation.gene_mu, perturbation.gene_sd);
        }
        else ds[i] = 1.0;
    }

    /* perturb transcript expression */
    double z = 0.0;
    for (i = 0; i < n; ++i) {
        phi[i] *= ds[gene_idx[i]];
        if (gsl_rng_uniform(rng) < perturbation.trans_pr) {
            phi[i] *= gsl_ran_lognormal(
                rng, perturbation.trans_mu, perturbation.trans_sd);
        }

        z += phi[i];
    }


    /* re-normalize */
    for (i = 0; i < n; ++i) phi[i] /= z;

    delete [] ds;
}


int seqsim_perturb(int argc, char* argv[])
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
                seqsim_perturb_help(stdout);                
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

    if (argc - opt_idx < 3) {
        seqsim_perturb_usage(stderr);
        return EXIT_FAILURE;
    }

    params P;
    P.read(argv[opt_idx++]);

    string perturb_name = argv[opt_idx++];
    if (P.perturbations.find(perturb_name) == P.perturbations.end()) {
        fprintf(stderr,
            "seqsim perturb: no perturbation named \"%s\" is given by the model.\n",
            perturb_name.c_str());
        return EXIT_FAILURE;
    }
    perturb_params& perturbation = P.perturbations[perturb_name];


    FILE* f = fopen(argv[opt_idx], "rb");
    if (f == NULL) {
        fprintf(stderr,
            "seqsim perturb: no expression file named \"%s\".\n", argv[opt_idx]);
        return EXIT_FAILURE;
    }

    vector<string> gene_ids;
    vector<string> trans_ids;
    vector<double> phi;

    map<string, double> expr;
    char  line[512];
    char *sep1, *sep2;
    while (fgets(line, sizeof(line), f)) {
        if ((sep1 = strchr(line, '\t')) == NULL) continue;
        if ((sep2 = strchr(sep1 + 1, '\t')) == NULL) continue;

        gene_ids.push_back(string(line, sep1 - line));
        trans_ids.push_back(string(sep1 + 1, sep2 - sep1 - 1));
        phi.push_back(atof(sep2 + 1));
    }

    fclose(f);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, rng_seed);
    
    perturb_phi(rng, perturbation, gene_ids, phi);

    gsl_rng_free(rng);

    size_t i;
    size_t n = phi.size();
    for (i = 0; i < n; ++i) {
        printf("%s\t%s\t%e\n",
            gene_ids[i].c_str(),
            trans_ids[i].c_str(),
            phi[i]);
    }

    return EXIT_SUCCESS;
}


