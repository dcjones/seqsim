/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef SEQSIM_FIT
#define SEQSIM_FIT

 #include <stdlib.h>

/* Fit a Gaussian mixture model.
 *    xs:  observations
 *    n:   number of observations
 *    k:   number of components
 *    mus: (output) k mean parameters
 *    sds: (output) k standard deviation parameters
 */
void fit_gaussian_mixture(
    const double* xs, size_t n, size_t k,
    double* ps, double* mus, double* sds);

#endif

