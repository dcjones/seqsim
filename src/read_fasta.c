
#include <string.h>
#include <stdbool.h>

#include "common.h"
#include "read_fasta.h"


/* Initialize an empty string. */
void fasta_str_init(fasta_str_t* str)
{
    str->n    = 0;
    str->size = 0;
    str->s    = NULL;
}


/* Free a string. */
void fasta_str_free(fasta_str_t* str)
{
    if (str) free(str->s);
}


/* Reserve space for `size` more characters. */
void fasta_str_reserve_extra(fasta_str_t* str, size_t size)
{
    if (str->n + size > str->size) {
        if (str->n + size > 2 * str->size) {
            str->size = str->n + size;
        }
        else str->size *= 2;
        str->s = realloc_or_die(str->s, str->size * sizeof(char));
    }
}


/* Append a char to a string. */
void fasta_str_append_char(fasta_str_t* a, char c)
{
    fasta_str_reserve_extra(a, 1);
    a->s[a->n++] = c;
}


/* Append some bytes to a string. */
void fasta_str_append(fasta_str_t* a, char* c, size_t n)
{
    fasta_str_reserve_extra(a, n);
    memcpy(a->s + a->n, c, n);
    a->n += n;
}


/* Read the next sequence in an open FASTA files.
 *
 * Args:
 *   file: An open file with a position before the next entry to be read.
 *   seq: A pointer to pointer which will be set to the sequence.
 *
 * Returns:
 *   true if a new sequence was read.
 *
 */
size_t read_fasta(FILE* file, namedseq_t** out)
{
    size_t size = 16;
    size_t n = 0;
    namedseq_t* seqs = malloc(size * sizeof(namedseq_t));

    const size_t bufsize = 1000000;
    char* buf = malloc_or_die(bufsize); buf[0] = '\0';
    char* next;
    size_t readlen;

    /* The three states are:
     *   0: Beginning of file.
     *   1. Reading a sequence name.
     *   2. Reading a sequence.
     */
    int state = 0;

    /* Are we at the beginning of a line. */
    bool linestart = true;

    while ((readlen = fread(buf, 1, bufsize, file))) {
        next = buf;
        char* end = buf + readlen;
        while (next < end) {
            /* found a new sequence */
            if (linestart && next[0] == '>') {
                state = 1;
                linestart = false;
                next++;

                if (n >= size) {
                    size *= 2;
                    seqs = realloc_or_die(seqs, size * sizeof(namedseq_t));
                }
                fasta_str_init(&seqs[n].name);
                fasta_str_init(&seqs[n].seq);
                ++n;

                continue;
            }

            char* u = memchr(next, '\n', end - next);
            if (u == NULL) {
                linestart = false;
                u = end;
            }
            else linestart = true;

            if (state == 1) {
                fasta_str_append(&seqs[n-1].name, next, u - next);
                if (linestart) {
                    fasta_str_append_char(&seqs[n-1].name, '\0');

                    /* ignore everything from the last ' '. */
                    char* c = strchr(seqs[n-1].name.s, ' ');
                    if (c) {
                        *c = '\0';
                        seqs[n-1].name.n = c - seqs[n-1].name.s;
                    }

                    fprintf(stderr, "Reading %s...\n", seqs[n-1].name.s);
                    state = 2;
                }
            }
            else if (state == 2) {
                fasta_str_append(&seqs[n-1].seq, next, u - next);
            }

            next = u + 1;
            linestart = true;
        }
    }

    free(buf);
    if(out) *out = seqs;
    return n;
}
