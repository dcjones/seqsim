
#ifndef SEQSIM_READ_FASTA

#if defined(__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* Simple encapsulated string. */
typedef struct
{
    char*          s;    /* null-terminated string */
    uint32_t       n;    /* length of s */
    uint32_t       size; /* bytes allocated for s */
} fasta_str_t;

void fasta_str_init(fasta_str_t* str);
void fasta_str_free(fasta_str_t* str);


/* A sequence with a name. */
typedef struct
{
    fasta_str_t name;
    fasta_str_t seq;
} namedseq_t;


size_t read_fasta(FILE* file, namedseq_t** out);

#if defined(__cplusplus)
}
#endif

#endif

