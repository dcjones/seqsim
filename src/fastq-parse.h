/*
 * This file is part of fastq-tools.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * parse :
 * A parser for FASTQ files.
 *
 */

#ifndef FASTQ_TOOLS_PARSE_H
#define FASTQ_TOOLS_PARSE_H

#if defined(__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <zlib.h>

/* share the gtf_parse.h definition of str_t */
#include "gtf_parse.h"

#if 0
typedef struct
{
    char*  s;    /* null-terminated string */
    size_t n;    /* length of s */
    size_t size; /* bytes allocated for s */
} str_t;
#endif



typedef struct
{
    str_t id1;
    str_t seq;
    str_t id2;
    str_t qual;
} seq_t;


seq_t* fastq_alloc_seq();
void fastq_free_seq(seq_t*);


typedef struct
{
    gzFile file;
    int    state;
    char*  buf;
    char*  c;
} fastq_t;


fastq_t* fastq_open(FILE*);
void fastq_close(fastq_t*);
int  fastq_next(fastq_t*, seq_t*);

void fastq_print(FILE* fout, seq_t* seq);


#if defined(__cplusplus)
}
#endif

#endif

