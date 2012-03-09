
#include "log_norm_mix.hpp"
#include <cmath>
#include "gsl/gsl_randist.h"

double log_norm_mix(
    gsl_rng* rng, size_t k, double* ps, double* mus, double* sds)
{
    double r = gsl_rng_uniform(rng);
    size_t i = 0;
    while (i < k && r > ps[i]) {
        r -= ps[i];
        ++i;
    }

    return exp(mus[i] + gsl_ran_gaussian(rng, sds[i]));
}

