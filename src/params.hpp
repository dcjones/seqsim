/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/* params:
 * Encapsulation of a set a parameters.
 */

#ifndef SEQSIM_PARAMS_HPP
#define SEQSIM_PARAMS_HPP

#include <stdint.h>
#include <stdlib.h>


struct params
{
    params();
    ~params();

    /* Read/write parameters from/to a YAML file. */
    void read(const char* fn);
    void write(const char* fn);

    /* Number of reads to generate. */
    uint64_t N;

    /* Generate paired-end reads. */
    bool paired;

    /* Strand-specificity: the probability the read is generated
     * in the same orientation as the transcript. */
    double strand_specificity;

    /* Log-normal mixture model from which gene expression is generated. */
    unsigned int gene_exp_k; /* number of components */
    double* gene_exp_p;
    double* gene_exp_mu;
    double* gene_exp_sd;

    /* Log-normal mixture model from which transcript expression is generated. */
    unsigned int trans_exp_k; /* number of components. */
    double* trans_exp_p;
    double* trans_exp_mu;
    double* trans_exp_sd;

    /* Size selection: fragments within the size bounds are generated.
     * If size_var is non-zero, these boundaries are fuzzy. Fragments
     * near the boundary are randomly accepted or rejected. */
    double size_lower;
    double size_upper;
    double size_var;
};


#endif

