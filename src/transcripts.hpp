/*
 * seqsim : the simplistic rna-seq simulator
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef SEQSIM_TRANSCRIPTS_HPP
#define SEQSIM_TRANSCRIPTS_HPP

#include "common.h"
#include "gtf_parse.h"
#include <string>
#include <vector>
#include <set>

struct exon
{
    exon(pos_t start, pos_t end);
    bool operator < (const exon& other) const;

    pos_t start;
    pos_t end;
};


struct transcript
{
    transcript();
    pos_t exonic_length() const;
    void add_exon(const gtf_row_t* row);

    strand_t strand;
    std::string seqname;
    std::string gene_id;
    std::string transcript_id;

    std::set<exon> exons;
};

void read_transcripts_from_gtf(FILE*, std::vector<transcript>& T);

#endif
