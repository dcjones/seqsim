#!/usr/bin/env python

# This is the first version of seqsim. I rewrote everything in C++ after, as the
# python version is just too slow.

'''
seqsim
------

This is an extremely simple RNA-Seq simulation. I make no guarantees of accuracy
or correctness. This is not a substitute for testing with real data!

Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
'''


# The 'gtf' module is my own. You need to get it from:
#        https://github.com/dcjones/cbgb
import gtf

import argparse
import Bio.SeqIO as SeqIO
import numpy as np
from numpy.random import \
    random, random_integers, binomial, geometric, \
    dirichlet, normal, shuffle

from sys         import stdout, stderr, stdin
from collections import defaultdict
from bisect      import bisect



def get_fasta_seqnames(fn):
    stderr.write('examining reference sequence ... ')

    seqnames = set()
    for line in open(fn):
        if line[0] == '>':
            seqnames.add(line[1:].strip('\n '))

    stderr.write('done. ({0} sequence found)\n'.format(len(seqnames)))
    return seqnames



def overlaps((u1, v1), (u2, v2)):
    ''' Determine whether [u1, v2] overlaps [u2, v2]. '''

    assert u1 <= v1
    assert u2 <= v2

    return u1 <= v2 and v1 >= v2


class transcript:
    ''' Representation of a single transcript. '''

    def __init__(self):
        self.exons = []
        self.seqname = None
        self.strand  = None

        self.gene_id       = None
        self.transcript_id = None


    def exonic_length(self):
        ''' Compute the length of the mature mRNA. '''

        return sum(end - start + 1 for (start, end) in self.exons)


    def add_exon(self, row):
        ''' Insert an exon into the transcript. '''

        if row.feature != 'exon': return

        if self.seqname is None: self.seqname = row.seqname
        else: assert self.seqname == row.seqname

        if self.strand is None: self. strand = row.strand
        else: assert self.strand == row.strand

        if self.gene_id is None: self.gene_id = row.attributes['gene_id']
        else: assert self.gene_id == row.attributes['gene_id']

        if self.transcript_id is None: self.transcript_id = row.attributes['transcript_id']
        else: assert self.transcript_id == row.attributes['transcript_id']


        # '-1' to make 0-based, end-inclusive
        (start, end) = (row.start - 1, row.end - 1)

        i = bisect(self.exons, (start, end))

        if len(self.exons) > i:
            assert not overlaps(self.exons[i], (start, end))

        if i > 0:
            assert not overlaps(self.exons[i - 1], (start, end))

        self.exons.insert(i, (start, end))




def read_transcripts(genes_fn, seqnames):
    ''' Parse all the transcripts from a GTF file. '''

    stderr.write('parsing gtf ... ')

    genes = defaultdict(transcript)

    for row in gtf.gtf_file(genes_fn):
        if row.seqname not in seqnames: continue

        genes[row.attributes['transcript_id']].add_exon(row)

    stderr.write('done. ({0} transcripts)\n'.format(len(genes)))
    return genes.values()



def print_expression(fout, T, xs):
    for (i, t) in enumerate(T):
        fout.write('{0}\t{1}\t{2}\n'.format(
                    t.gene_id, t.transcript_id, xs[i]))



def generate_expression(n):
    '''
    Generate random expression data in terms of proportion of RNA molecules.
    '''

    stderr.write('generating expression ... ')

    # choose the number of expressed transcripts
    m = binomial(n, args.expr_pr)

    # choose which transcripts are expressed
    us = np.array(range(n))
    shuffle(us)
    us = us[:m]

    # generate expression values
    ys = dirichlet(np.repeat(args.expr_a, m))
    xs = np.zeros(n)
    for (i, j) in enumerate(us):
        xs[j] = ys[i]

    stderr.write('done.\n')

    return xs


def generate_fragments(T, xs, m):
    '''
    Generate m random fragments from the transcripts T.

     The model of random fragmentation used here is very simple: we simple
     choose two random break-points, as follows,

     1. Choose a transcript t randomly, with probability proportional to its
        relative abundance times its length.

     2. Choose a start breakpoint uniformly from within the mature mRNA.

     3. Choose an end breakpoint by drawing a geometric variate and adding it
        to the start breakpoint. If this ends up outside the mRNA, reject, and
        try again from step 1.

    '''


    stderr.write('generating fragments ...\n')

    ws = np.array([t.exonic_length() for t in T]) * xs
    ws /= sum(ws)
    np.cumsum(ws, out = ws)

    fs = np.empty((m, 3), dtype = int)

    i = 0
    while i < m:
        # choose a random transcript
        j = ws.searchsorted(random())
        l = T[j].exonic_length()

        # choose start and end points
        start = random_integers(0, l - 1)
        end   = start + geometric(args.frag_pr)

        # reject improper fragments (end breakpoint outside the mRNA)
        if end >= l: continue

        # reject fragments lost in size selection
        fraglen = end - start + 1
        if fraglen < args.size_low + normal(args.size_noise):
            continue

        if fraglen > args.size_high + normal(args.size_noise):
            continue

        # reject fragments too short to sequence
        if fraglen < args.readlen:
            continue


        fs[i] = (j, start, end)
        i += 1

        if (i % 100000) == 0:
            stderr.write('\t{0}\n'.format(i))

    stderr.write('done. ({0} fragments)\n'.format(m))

    return fs


def generate_reads(fs, m):
    '''
    Given size-selected fragments, generated m reads.

    This is done by simple sampling fragments (with replacement), and returning
    the indexes.
    '''

    stderr.write('generating reads ... ')

    rs = np.zeros(len(fs), dtype = int)

    for i in xrange(m):
        i = random_integers(0, len(fs) - 1)
        rs[i] += 1

    stderr.write('done.\n')

    return rs




def get_seq(start, end, t, seq):
    '''
    Get the genomic sequence for positions 'start'
    through 'end from within transcript t.
    '''
    s = None
    z = 0

    gene_start = t.exons[0][0]

    for exon in t.exons:
        l = exon[1] - exon[0] + 1

        if z <= start < z + l:
            ss = seq[gene_start + max(start, z) : \
                     gene_start + min(end + 1, z + l)]

            if s is None: s = ss
            else:         s += ss

            start = z + l
            if start > end: break

        z += l

    return s



def print_read(outf, readnum, start, end, t, seq):
    '''
    Extract the sequence and print a read in FASTQ format.
    '''

    seq1 = get_seq(start, start + args.readlen - 1, t, seq)
    seq2 = get_seq(end - args.readlen + 1, end, t, seq)


    if args.stranded:
        strand = t.strand
    else:
        if random() < 0.5:
            strand = '+'
        else:
            strand = '-'


    if strand == '+':
        seq2 = seq2.reverse_complement()
    else:
        (seq1, seq2) = (seq2.reverse_complement(), seq1)


    assert len(seq1) == args.readlen
    assert len(seq2) == args.readlen

    read_id = 'seqsim.{prefix}.{readnum:09d}'.format(
            prefix  = args.out,
            readnum = readnum)

    if args.paired_end:
        outf.write('@{read_id}\n{seq}\n+\n{qual}\n'.format(
            read_id = read_id + '/1',
            seq     = seq1.tostring().upper(),
            qual    = 'h' * args.readlen))

        outf.write('@{read_id}\n{seq}\n+\n{qual}\n'.format(
            read_id = read_id + '/2',
            seq     = seq1.tostring().upper(),
            qual    = 'h' * args.readlen))
    else:
        outf.write('@{read_id}\n{seq}\n+\n{qual}\n'.format(
            read_id = read_id,
            seq     = seq1.tostring().upper(),
            qual    = 'h' * args.readlen))



def main():
    ap = argparse.ArgumentParser()

    ap.add_argument('-n', type = int, default = 10000000,
                    help = 'number of reads / pairs to generate')

    ap.add_argument('-o', '--out', type = str, default = 'seqsim',
                    help = 'prefix for output files')

    ap.add_argument('-p', '--paired-end', action = 'store_true',
                    default = False,
                    help = 'generate paired end reads')

    ap.add_argument('-s', '--stranded', action = 'store_true',
                    default = False,
                    help = 'data generated in strand specific')

    ap.add_argument('-k', '--readlen', type = int, default = 75,
                    help = 'length of each read to generate')

    ap.add_argument('--expr_pr', type = float, default = 0.4,
                    help = 'the probability a transcript is expressed')

    ap.add_argument('--expr_a', type = float, default = 1.0,
                    help = 'parameter to the symmetricx dirichlet distribution'
                           'from which expression value are drawn')

    ap.add_argument('--frag_pr', type = float, default = 0.004,
                    help = 'fragmentation probability')

    ap.add_argument('--frag_n', type = int, default = 100000000,
                    help = 'number of fragments to generate')

    ap.add_argument('--size_low', type = float, default = 180,
                    help = 'low threshold for size selection')

    ap.add_argument('--size_high', type = float, default = 220,
                    help = 'low threshold for size selection')

    ap.add_argument('--size_noise', type = float, default = 220,
                    help = 'parameter controlling how inexact size selection is')

    ap.add_argument('genome_fn', metavar = 'genome.fa', type = str,
                    help = 'genome sequence, in FASTA format')

    ap.add_argument('genes_fn', metavar = 'genes.gtf', type = str,
                    help = 'gene annotations in GTF format')

    global args
    args = ap.parse_args()

    seqnames = get_fasta_seqnames(args.genome_fn)

    T = read_transcripts(args.genes_fn, seqnames)
    n = len(T)

    # generate / output true expression
    xs = generate_expression(len(T))
    print_expression(open(args.out + '.expr', 'w'), T, xs)

    # generate random fragments
    fs = generate_fragments(T, xs, args.frag_n)

    # sample reads from random fragments
    rs = generate_reads(fs, args.n)

    # extract sequence
    stderr.write('getting sequence ...\n')
    fout = open(args.out + '.fastq', 'w')

    readnum = 1
    for seq in SeqIO.parse(open(args.genome_fn), 'fasta'):
        stderr.write('\t{seqname}\n'.format(seqname = seq.name))

        # indicies of transcripts from this sequence
        us = set(u for u in xrange(n) if T[u].seqname == seq.name)

        # indices of fragments from this sequence
        vs = [v for v in xrange(args.frag_n) if fs[v, 0] in us]

        for v in vs:
            for _ in xrange(rs[v]):
                print_read(fout, readnum, fs[v, 1], fs[v, 2], T[fs[v, 0]], seq.seq)
                readnum += 1

    stderr.write('done.\n')



if __name__ == '__main__': main()


