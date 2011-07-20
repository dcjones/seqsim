
#ifndef ISOLATOR_COMMON_H
#define ISOLATOR_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

/* genomic position */
typedef long pos_t;

/* sequence identifier */
typedef int seqid_t;

typedef enum {
    strand_pos = 0,
    strand_neg = 1,
    strand_na  = 2
}strand_t;

void or_die(int b, const char* msg);

void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);
FILE* fopen_or_die(const char* fn, const char* mode);


#ifdef __cplusplus
}
#endif

#endif

