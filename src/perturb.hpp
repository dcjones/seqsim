/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef SEQSIM_PERTURB_HPP
#define SEQSIM_PERTURB_HPP

#include <cstdio>

void seqsim_perturb_usage(FILE*);
void seqsim_perturb_help(FILE*);

int seqsim_perturb(int argc, char* argv[]);

#endif
