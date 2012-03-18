
/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "config.h"
#include "params.hpp"
#include "express.hpp"
#include "perturb.hpp"
#include "generate.hpp"
#include "common.h"
#include "gtf_parse.h"
#include "transcripts.hpp"
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


static void print_version(FILE* fout)
{
    fprintf(fout, "seqsim %s\n", VERSION);
}


static void print_usage(FILE* fout)
{
    fprintf(fout,
        "Usage: seqsim <command> [<arguments>]\n\n"
        "Where <command> is one of:\n"
        "    express\n"
        "    perturb\n"
        "    generate\n\n"
        "See 'seqsim help <command>' for more.\n"
        );
}


int main(int argc, char* argv[])
{
    if (argc < 2) {
        print_usage(stdout);
        return EXIT_SUCCESS;
    }

    argv++;
    argc--;

    // special case for --version
    if (strcmp(argv[0], "--version") == 0) {
        print_version(stdout);
        return EXIT_SUCCESS;
    }
    else if (strcmp(argv[0], "express") == 0) {
        return seqsim_express(argc, argv);
    }
    else if (strcmp(argv[0], "perturb") == 0) {
        return seqsim_perturb(argc, argv);
    }
    else if (strcmp(argv[0], "generate") == 0) {
        return seqsim_generate(argc, argv);
    }

    fprintf(stderr, "Unknown command %s.\n\n", argv[0]);
    return EXIT_FAILURE;
}


