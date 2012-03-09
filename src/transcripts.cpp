
#include "transcripts.hpp"
#include <map>

using namespace std;

exon::exon(pos_t start, pos_t end)
    : start(start)
    , end(end)
{}

bool exon::operator < (const exon& other) const
{
    if (start != other.start) return start < other.start;
    else                      return end   < other.end;
}


transcript::transcript()
    : strand(strand_na)
{
}


pos_t transcript::exonic_length() const
{
    pos_t len = 0;
    set<exon>::iterator i;
    for (i = exons.begin(); i != exons.end(); ++i) {
        len += i->end - i->start + 1;
    }

    return len;
}


void transcript::add_exon(const gtf_row_t* row)
{
    if (strcmp(row->feature->s, "exon") != 0) return;

    str_t* val;

    if (seqname.empty())     seqname = row->seqname->s;
    if (strand == strand_na) strand  = row->strand;
    if (gene_id.empty()) {
        val = (str_t*) str_map_get(row->attributes, "gene_id", 7);
        gene_id = val->s;
    }

    if (transcript_id.empty()) {
        val = (str_t*) str_map_get(row->attributes, "transcript_id", 13);
        transcript_id = val->s;
    }

    exons.insert(exon(row->start - 1, row->end - 1));
}


void read_transcripts_from_gtf(FILE* fin, std::vector<transcript>& T)
{
    fprintf(stderr, "parsing gtf ... ");

    map<string, transcript> ts;

    gtf_file_t* gf = gtf_file_alloc(fin);
    gtf_row_t* row = gtf_row_alloc();
    str_t* val;

    map<string, transcript>::iterator i;
    while (gtf_next(gf, row)) {
        val = (str_t*) str_map_get(row->attributes, "transcript_id", 13);
        if (val && val->s) ts[string(val->s)].add_exon(row);
    }

    gtf_row_free(row);
    gtf_file_free(gf);

    T.reserve(ts.size());
    for (i = ts.begin(); i != ts.end(); ++i) {
        T.push_back(i->second);
    }

    fprintf(stderr, "done.\n");
}

