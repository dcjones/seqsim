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
#include <cstdlib>
#include <map>
#include <vector>
#include <string>


struct perturb_params
{
    perturb_params();
    ~perturb_params();

    // Probability that a gene is effected by the perturbation.
    double gene_pr;

    // Normal distribution parameters by which gene expression is permuted
    double gene_mu;
    double gene_sd;

    // Probability that a transcript is effected by the perturbation.
    double trans_pr;

    // Normal distribution parameters by which gene expression is permuted
    double trans_mu;
    double trans_sd;
};


struct params
{
    params();
    ~params();

    /* Read/write parameters from/to a YAML file. */
    void read(const char* fn);
    void write(const char* fn);

    /* Number of reads to generate. */
    uint64_t N;

    /* read length */
    uint32_t readlen;

    /* Generate paired-end reads. */
    bool paired;

    /* Strand-specificity: the probability the read is generated
     * in the same orientation as the transcript. */
    double strand_specificity;

    /* Log-normal mixture model from which gene expression is generated. */
    unsigned int gene_exp_k; /* number of components */
    std::vector<double> gene_exp_p;
    std::vector<double> gene_exp_mu;
    std::vector<double> gene_exp_sd;

    /* Symmetric Dirichlet distribution from which transcript expression is generated. */
    double trans_exp_alpha;

    /* Size selection: fragments within the size bounds are generated.
     * If size_var is non-zero, these boundaries are fuzzy. Fragments
     * near the boundary are randomly accepted or rejected. */
    double size_lower;
    double size_upper;
    double size_sd;

    /* Expression perturbations */
    std::map<std::string, perturb_params> perturbations;
};


#endif

