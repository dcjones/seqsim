
#include "generate.hpp"
#include <cstdlib>
#include <getopt.h>

void seqsim_generate_usage(FILE* fout)
{
    fprintf(fout,
        "Usage: seqsim generate [options] <genome.fa> <genes.gtf> <expression.yaml>\n");
}


void seqsim_generate_help(FILE* fout)
{
    seqsim_generate_usage(fout);
    fprintf(fout, "TODO");
}


int seqsim_generate(int argc, char* argv[])
{
    /* TODO */

    return EXIT_SUCCESS;
}

