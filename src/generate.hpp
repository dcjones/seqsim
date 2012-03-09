/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef SEQSIM_GENERATE_HPP
#define SEQSIM_GENERATE_HPP

#include <cstdio>

void seqsim_generate_usage(FILE*);
void seqsim_generate_help(FILE*);

int seqsim_generate(int argc, char* argv[]);

#endif
