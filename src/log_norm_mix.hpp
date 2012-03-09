/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef SEQSIM_LOG_NORM_MIX_HPP
#define SEQSIM_LOG_NORM_MIX_HPP

#include <gsl/gsl_rng.h>

/* Draw a random value from a log-normal mixture model, where
 * ps are the mixing coefficients, mus are the means, and sds
 * are the standard deviations.
 */
double log_norm_mix(gsl_rng* rng, size_t k, double* ps, double* mu, double* sd);

#endif


