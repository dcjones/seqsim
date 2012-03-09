/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef SEQSIM_EXPRESS_HPP
#define SEQSIM_EXPRESS_HPP

#include <cstdio>

void seqsim_express_usage(FILE*);
void seqsim_express_help(FILE*);

int seqsim_express(int argc, char* argv[]);

#endif
