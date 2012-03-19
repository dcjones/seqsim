

#include "gaussian_mixture.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int main(int argc, char* argv[])
{
    FILE* fin;

    if (argc > 1) {
        fin = fopen(argv[1], "r");
        if (fin == NULL) {
            fprintf(stderr, "fit_lnorm_mix: can't open file \"%s\".\n", argv[1]);
            return EXIT_FAILURE;
        }
    }
    else fin = stdin;


    size_t k = 2; /* TODO: option to set this */
    size_t n = 0;
    size_t size = 10000;
    char* endptr;
    double x;
    double* xs = malloc(size * sizeof(double));

    char line[512];
    while (fgets(line, sizeof(line), fin)) {
        if (n == size) {
            size += 1000;
            xs = realloc(xs, size * sizeof(double));
        }

        x = strtod(line, &endptr);
        if (x == 0.0 || endptr == line) continue;

        xs[n++] = log(x);
    }

    if (fin != stdin) fclose(fin);

    fprintf(stderr, "read %zu values.\n", n);
    fprintf(stderr, "fitting ... \n");

    double* ps  = malloc(k * sizeof(double));
    double* mus = malloc(k * sizeof(double));
    double* sds = malloc(k * sizeof(double));

    fit_gaussian_mixture(xs, n, k, ps, mus, sds);
    
    free(ps);
    free(mus);
    free(sds);
    free(xs);

    return EXIT_SUCCESS;
}




