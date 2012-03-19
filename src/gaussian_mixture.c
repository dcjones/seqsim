#include "gaussian_mixture.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

static const double sqrt_2pi = 2.5066282746310002;
static double sq(double x) { return x * x; }
static double pnorm(double x, double mu, double sd)
{
    return exp(-sq(x - mu) / (2.0 * sq(sd))) / (sqrt_2pi * sd);
}

static void e_step(
    size_t n, size_t k,
    const double* xs, double* zs,
    const double* ps, const double* mus, const double* sds)
{
    size_t i, j;
    for (i = 0; i < k; ++i) {
        for (j = 0; j < n; ++j) {
            zs[i * n + j] = ps[i] * pnorm(xs[j], mus[i], sds[i]);
        }
    }

    double z_sum;
    for (j = 0; j < n; ++j) {
        z_sum = 0.0;
        for (i = 0; i < k; ++i) z_sum += zs[i * n + j];
        for (i = 0; i < k; ++i) zs[i * n + j] /= z_sum;
    }
}


static void m_step(
    size_t n, size_t k,
    const double* xs, const double* zs,
    double* ps, double* mus, double* sds)
{
    size_t i, j;
    double z_sum;
    for (i = 0; i < k; ++i) {
        z_sum = 0.0;
        mus[i] = 0.0;
        for (j = 0; j < n; ++j) {
            z_sum += zs[i * n + j];
            mus[i] += zs[i * n + j] * xs[j];
        }
        mus[i] /= z_sum;

        sds[i] = 0.0;
        for (j = 0; j < n; ++j) {
            sds[i] += zs[i * n + j] * sq(xs[j] - mus[i]);
        }
        sds[i] = sqrt(sds[i] / z_sum);

        ps[i] = z_sum / (double) n;
    }
}


static void print_progress(
        size_t k,
        const double* ps,
        const double* mus,
        const double* sds)
{
    size_t i;

    fprintf(stderr, "------------------------------------\n");

    fprintf(stderr, "ps  = [ %e", ps[0]);
    for (i = 1; i < k; ++i) fprintf(stderr, ", %e", ps[i]);
    fprintf(stderr, " ]\n");

    fprintf(stderr, "mus = [ %e", mus[0]);
    for (i = 1; i < k; ++i) fprintf(stderr, ", %e", mus[i]);
    fprintf(stderr, " ]\n");

    fprintf(stderr, "sds = [ %e", sds[0]);
    for (i = 1; i < k; ++i) fprintf(stderr, ", %e", sds[i]);
    fprintf(stderr, " ]\n");
}


void fit_gaussian_mixture(
    const double* xs, size_t n, size_t k,
    double* ps, double* mus, double* sds)
{
    /* expectations of latent indicator variables assigning
     * observations to mixture components */
    double* zs = malloc(k * n * sizeof(double));

    /* set initial mixture probabilities to uniform */
    size_t i;
    double init_p = 1.0 / (double) k;
    for (i = 0; i < k; ++i) ps[i] = init_p;

    /* set initial component std. deviations to the
     * sample std. deviation. */
    double mu = 0.0;
    for (i = 0; i < n; ++i) mu += xs[i];
    mu /= (double) n;

    double sd = 0.0;
    for (i = 0; i < n; ++i) sd += sq(xs[i] - mu);
    sd /= (double) (n - 1);

    for (i = 0; i < k; ++i) sds[i] = sd;

    /* set initial evenly spaced component means */
    for (i = 0; i < k; ++i) {
        mus[i] = mu + 3.0 * sd * (2.0 * ((double) i / (double) k) - 1.0);
    } 

    const double eps = 1e-8;
    bool converged;

    double* mus0 = malloc(k * sizeof(double));
    double* sds0 = malloc(k * sizeof(double));

    do {
        print_progress(k, ps, mus, sds);
        memcpy(mus0, mus, k * sizeof(double));
        memcpy(sds0, sds, k * sizeof(double));

        e_step(n, k, xs, zs, ps, mus, sds);
        m_step(n, k, xs, zs, ps, mus, sds);

        converged = true;
        for (i = 0; i < k; ++i) {
            if (fabs(mus[i] - mus0[i]) > eps || fabs(sds[i] - sds0[i]) > eps) {
                converged = false;
                break;
            }
        }

    } while(!converged);

    free(zs);
    free(mus0);
    free(sds0);
}



