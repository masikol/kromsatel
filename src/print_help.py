# -*- coding: utf-8 -*-

import sys

def print_help(version, last_update_data):
    print('kromsatel')
    print('Version {}. {} edition.'.format(version, last_update_data))
    print("\nUsage:")
    print('  ./kromsatel.py <YOUR_READS> -d <DATABASE_WITH_FRAGMENTS>')
    print('\nOptions:')
    print('  -h (--help) -- print help message and exit.')
    print('  -v (--version) -- print version and exit.')
    print("""  -o (--outdir) -- output directory.
    Default value is `./kromsatel_output`""")
    print("""  -d (--db) -- path to BLAST database of amplicons.
    Mandatory option (it sounds like an oxymoron but anyway).""")
    print("""  -p (--primers-to-rm) -- CSV file containing primers and their names.
    It must be specified if you intend to remove primer sequences from reads.
    See section "Removing primer sequences" for details.""")
    print("""  -t (--threads) -- number of threads to launch.
    Default: 1 thread.""")
    print("""  --am -- minimum length of alignment against major amplicon.
    Shorter alignments will not be considered.
    Default: 100 bp.""")
    print("""  --im -- minimum length of alignment against minor amplicon.
    Shorter alignments will not be considered.
    Default: 25 bp.""")
    print("""  -c (--chunk-size) -- number of reads to blast within a single query.
    Default: 1000 reads.""")
    print('\nExample:')
    print('  ./kromsatel.py corona_reads.fastq -d fragments-db/nCoV-2019_ref-fragments.fasta')
# end def print_help
