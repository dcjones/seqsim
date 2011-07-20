
#include "common.h"
#include <stdio.h>


void or_die(int b, const char* msg)
{
    if (b == 0) {
        fputs(msg, stderr);
        exit(1);
    }
}


void* malloc_or_die(size_t n)
{
    void* ptr = malloc(n);
    if (ptr == NULL) {
        fprintf(stderr, "Falied to allocate %zu bytes. Out of memory.\n", n);
        exit(EXIT_FAILURE);
    }

    return ptr;
}




void* realloc_or_die(void* ptr, size_t n)
{
    ptr = realloc(ptr, n);
    if (ptr == NULL) {
        fprintf(stderr,"Falied to (re)allocate %zu bytes. Out of memory.\n", n);
        exit(EXIT_FAILURE);
    }

    return ptr;
}



FILE* fopen_or_die(const char* fn, const char* mode)
{
    FILE* f = fopen(fn, mode);
    if (f == NULL) {
        fprintf(stderr, "Cannot open file %s.\n", fn);
        exit(EXIT_FAILURE);
    }

    return f;
}


