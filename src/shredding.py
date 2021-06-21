# -*- coding: utf-8 -*-

import io
import os
import sys
import tempfile
import operator
import functools
import multiprocessing as mp

import src.blast
import src.fastq
import src.filesystem
from src.filesystem import OPEN_FUNCS
from src.platform import platf_depend_exit
from src.printing import getwt

try:
    import numpy as np
except ImportError:
    print('Error: package `numpy` is not installed')
    print('Please, install it. Example command:')
    print('  pip3 install numpy')
    platf_depend_exit(1)
# end try

try:
    import pandas as pd
except ImportError:
    print('Error package `pandas` is not installed')
    print('Please, install it. Example command:')
    print('  pip3 install pandas')
    platf_depend_exit(1)
# end try


# Locks and "synchronized" values
WRITE_LOCK = mp.Lock()             # Lock for writing to output file
PRINT_LOCK = mp.Lock()             # Lock for printing status bar
INC_LOCK = mp.Lock()               # Lock for incrementing INC_VAL
INC_VAL = None                     # Counter for counting processed reads
                                   # (must be initialized just before pool start)
NEXT_PRINT_NUM = None              # Next number of processed reads to update status bar
                                   # (must be initialized just before pool start)


def _str2df(align_result_str):
    # Function converts tabular string data to pandas.DataFrame.
    #
    #:param align_result_str: tabular string to be converted to pandas.DataFrame;
    #:type align_result_str: str;
    #
    # Returns pandas.DataFrame contatining the same information as the string passed.

    # aln_df = pd.DataFrame([l.split('\t') for l in align_result_str.split('\n')])
    aln_df = pd.read_csv(io.StringIO(align_result_str),
        sep='\t',
        header=None, # blastn returns no header in "outfmt 6"
        names= ['qseqid', 'sseqid', 'sstrand',
                'length',   'qlen',    'slen',
                'qstart',   'qend',  'sstart',
                                       'send'],
        true_values=['plus'], # plus strand will be True
        false_values=['minus'] # minus strand will be False
    )

    return aln_df
# end def _str2df


def _rm_spuriously_major_alns(major_sorted_alns, minor_sorted_alns, primers_lengths):
    # Function removes spuriously "major" alignments from `major_sorted_alns`
    #   according to `minor_sorted_alns` and `primers_lengths` in order to
    #   remove primers sequences successfully.
    # :param major_sorted_alns: dataframe with "major" alignments;
    # :type major_sorted_alns: pandas.DataFrame;
    # :param minor_sorted_alns: dataframe with "minor" alignments;
    # :type minor_sorted_alns: pandas.DataFrame;
    # :param primers_lengths: tuple mapping primer's index to it's length.
    #   In this tuple, index of a primer equals index of this primer in source .csv file;
    # :type primers_lengths: tuple<int>;
    # Returns `major_sorted_alns` with removed spurious aligments.

    # Initizlize list for IDs (in DataFrame index) of major alignments to remove.
    major_idx_to_rm = list()

    # Iterate over major alignments
    for i_major, major_aln in major_sorted_alns.iterrows():

        # Find minor alignments starting at the same position (within query sequence)
        #   as current `major_aln`
        qstart = major_aln['qstart']
        minor_alns_with_same_start = minor_sorted_alns[minor_sorted_alns['qstart'] == qstart]

        # If there are any (there must be, actually), the first one must be ours
        if minor_alns_with_same_start.shape[0] != 0:

            # Store this minor alignment in a local variable
            minor_aln = minor_alns_with_same_start.iloc[0]

            # Get index of the amplicon against which `minor_aln` is aligned
            forw_primer_amplicon_idx = int(minor_aln['sseqid'][1:])

            # Obtain index of the reverse primer for this minor amplicon,
            #   considering `sstrand`
            reverse_primer_idx = int(2 * forw_primer_amplicon_idx - 1)
            if minor_aln['sstrand'] == False:
                reverse_primer_idx += 1
            # end if
            reverse_primer_len = primers_lengths[reverse_primer_idx]

            # Check if `major_aln` is long enough to be kept
            len_diff = abs(major_aln['send'] - major_aln['sstart']) \
                       - minor_aln['slen']

            # If it is no longer than `minor_aln['slen'] + reverse_primer_len`,
            #   then we'll remove it
            if len_diff <= reverse_primer_len:
                major_idx_to_rm.append(i_major)
                continue
            # end if
        # end if

        # Find minor alignments ending at the same position (within query sequence)
        #   as current `major_aln`
        qend = major_aln['qend']
        minor_alns_with_same_end = minor_sorted_alns[minor_sorted_alns['qend'] == qend]

        # If there are any (there must be, actually), the first one must be ours
        if minor_alns_with_same_end.shape[0] != 0:

            # Store this minor alignment in a local variable
            minor_aln = minor_alns_with_same_end.iloc[0]

            # Get index of the amplicon against which `minor_aln` is aligned
            forw_primer_amplicon_idx = int(minor_aln['sseqid'][1:])

            # Obtain index of the reverse primer (in notation of primers file, it'll be forward)
            #   for this minor amplicon, considering `sstrand`
            forward_primer_idx = int(2 * forw_primer_amplicon_idx)
            if minor_aln['sstrand'] == False:
                forward_primer_idx -= 1
            # end if
            forward_primer_len = primers_lengths[forward_primer_idx]

            # Check if `major_aln` is long enough to be kept
            len_diff = abs(major_aln['send'] - major_aln['sstart']) \
                       - minor_aln['slen']

            # If it is no longer than `minor_aln['slen'] + forward_primer_len`,
            #   then we'll remove it
            if len_diff <= forward_primer_len:
                major_idx_to_rm.append(i_major)
                continue
            # end if
        # end if
    # end for

    # Remove spurious "major" amplicons
    major_sorted_alns.drop(major_idx_to_rm, inplace=True)

    return major_sorted_alns
# end def _rm_spuriously_major_alns


def _get_aligned_spans_rm_primers(curr_alns, primers_lengths, min_len_minor):
    # Function analyses obtained alignments and find "spans" -- subsequences
    #   into which input read should be "shredded".
    # These spans should not overlap.
    # The function removes primers sequences, which might by left
    #   due to "spurious major" alignments. Also, the function removes
    #   short (`min_len_minor`) "minor" alignments.
    #
    # :param curr_alns: dataframe containing alignmens data obtained from _str2df function;
    # :type curr_alns: pandas.DataFrame;
    # :param primers_lengths: tuple mapping primer's index to it's length.
    #   In this tuple, index of a primer equals index of this primer in source .csv file;
    # :type primers_lengths: tuple<int>;
    # :param min_len_minor: minimum length of "minor" alignment;
    # :type min_len_minor: int;

    # If there are no alignments -- just return empty list.
    if curr_alns.empty:
        return tuple()
    # end if

    curr_alns_major_sorted = curr_alns[curr_alns['major'] == True]\
        .sort_values(by='length', ascending=False)
    curr_alns_minor_sorted = curr_alns[curr_alns['major'] == False]\
        .sort_values(by='length', ascending=False)

    try:
        # Remove spuriously "major" alignments
        curr_alns_major_sorted = _rm_spuriously_major_alns(
            curr_alns_major_sorted,
            curr_alns_minor_sorted,
            primers_lengths
        )
    except ValueError as err:
        print('Error: {}'.format(err))
        print("The reason is probably that you've changed titles of the amplicons \
in fasta file befor creating the database.")
        platf_depend_exit(1)
    # end try

    # List of result spans
    aligned_spans = list()

    # Select non-overlapping alignments and remove short "minor" alignments

    does_overlap = lambda span: (span[0] <= aln['qstart']-1 <= span[1])\
                             or (span[0] <= aln['qend']     <= span[1])

    for aln_collection in (curr_alns_major_sorted, curr_alns_minor_sorted):
        for _, aln in aln_collection.iterrows():

            # If we have short "minor" alignment -- omit it
            if aln['major'] == False and aln['length'] < min_len_minor:
                continue
            # end if

            overlap_result = functools.reduce(operator.or_,
                map(does_overlap, aligned_spans),
                False
            )

            if not overlap_result:
                aligned_spans.append( (aln['qstart']-1, aln['qend']) )
            # end if
        # end for

    return aligned_spans
# end def _get_aligned_spans_rm_primers


def _get_aligned_spans_keep_primers(curr_alns):
    # Function analyses obtained alignments and find "spans" -- subsequences
    #   into which input read should be "shredded".
    # These spans should not overlap.
    # The function pays no attention on primer sequences.
    #
    # :param curr_alns: dataframe containing alignmens data obtained from _str2df function;
    # :type curr_alns: pandas.DataFrame;

    # If there are no alignments -- just return empty list.
    if curr_alns.empty:
        return tuple()
    # end if

    curr_alns_major_sorted = curr_alns[curr_alns['major'] == True]\
        .sort_values(by='length', ascending=False)
    curr_alns_minor_sorted = curr_alns[curr_alns['major'] == False]\
        .sort_values(by='length', ascending=False)

    # List of result spans
    aligned_spans = list()

    # Select non-overlapping alignments

    does_overlap = lambda span: (span[0] <= aln['qstart']-1 <= span[1])\
                             or (span[0] <= aln['qend']     <= span[1])

    for aln_collection in (curr_alns_major_sorted, curr_alns_minor_sorted):
        for _, aln in aln_collection.iterrows():

            overlap_result = functools.reduce(operator.or_,
                map(does_overlap, aligned_spans),
                False
            )

            if not overlap_result:
                aligned_spans.append( (aln['qstart']-1, aln['qend']) )
            # end if
        # end for

    return aligned_spans
# end def _get_aligned_spans_keep_primers


def _get_bar_len():
    # Function returns length of status bar.
    # Firstly it tries to obtain "adaptive" bar length from terminal width.
    # It might fail if output is redirected somewhere, e.g. with `tee`.
    # If obtaining "adaptive" length fails, bar length will be set to some low value.
    # Return type: int.
    try:
        bar_len = int(os.get_terminal_size().columns * 0.40)
    except OSError:
        # Default low value
        bar_len = 40
    # end try
    return bar_len
# end _get_bar_len


def _set_major(row):
    # Function to be used with pandas.DataFrame.apply function.
    # It adds a column indicating if `sseqid` of the alignment
    #   in a particular row is major fragment, i.e. it starts with 'A'.
    #
    # :param row: row of dataframe to which this function if applied;
    # :type row: pandas.Series;

    if row['sseqid'][0] == 'A':
        row['major'] = True
    # end if

    return row
# end def _set_major


def _add_major_col(aln_df):
    # And a column indicating if subject sequence is a major fragment
    # :param aln_df: data frame of alignments;
    # :type aln_df: pandas.DataFrame;

    aln_df.insert(loc=aln_df.shape[1],
        column='major',
        value=np.repeat(False, aln_df.shape[0])
    )
    aln_df = aln_df.apply(_set_major, axis=1)

    return aln_df
# end def _add_major_col


def _filter_short_alns(aln_df, min_len_major, min_len_minor):
    # Function removes too short alignments from a data frame,
    #   cosidering if alignment is "major" or "minor"
    #
    # :param aln_df: data frame of alignments;
    # :type aln_df: pandas.DataFrame;
    # :param min_len_major: minimum length of "major" alignment;
    # :type min_len_major: int;
    # :param min_len_minor: minimum length of "minor" alignment;
    # :type min_len_minor: int;

    aln_df = aln_df[
        ((aln_df['major'] == True)  & (aln_df['length'] >= min_len_major))\
      | ((aln_df['major'] == False) & (aln_df['length'] >= min_len_minor))
    ]
    return aln_df
# end def _filter_short_alns


def _filter_short_major_alns(aln_df, min_len_major):
    # Function removes too short "major" alignments from a data frame.
    #
    # :param aln_df: data frame of alignments;
    # :type aln_df: pandas.DataFrame;
    # :param min_len_major: minimum length of "major" alignment;
    # :type min_len_major: int;

    aln_df = aln_df[
        (aln_df['major'] == False)
        | ((aln_df['major'] == True)  & (aln_df['length'] >= min_len_major))
    ]
    return aln_df
# end def _filter_short_major_alns


def _shredder(
    fq_chunk,
    db_fpath,
    min_len_major,
    min_len_minor,
    nreads,
    inc_num,
    outfpath,
    primers_lengths
    ):
    # "Kernel" function, which alignes and splits chimeric reads into non-chimeric fragments.
    # :param fq_chunk: current fastq chunk to process;
    # :type fq_chunk: dict<str: dict<str: str>>;
    # :param db_fpath: path to BLAST database to use;
    # :type db_fpath: str;
    # :param min_len_major: minimum length of "major" alignment;
    # :type min_len_major: int;
    # :param min_len_minor: minimum length of "minor" alignment;
    # :type min_len_minor: int;
    # :param nreads: total number of reads in current file;
    # :type nreads: int;
    # :param inc_num: increment value for status bar;
    # :type inc_num: int;
    # :param outfpath: path to output file;
    # :type outfpath: str;
    # :param primers_lengths: tuple mapping primer's index to it's length.
    #   In this tuple, index of a primer equals index of this primer in source .csv file;
    # :type primers_lengths: tuple<int>;

    # Obtain path to tempfile for writing auries into it
    query_fpath = os.path.join(
        tempfile.gettempdir(),
        'kromsatel_query_{}.fasta'.format(os.getpid())
    )

    # Configure command line for BLASTn
    blast_cmd = src.blast.configure_blastn_cmd(query_fpath, db_fpath)

    # Convert input fastq file to fasta format in order to pass the latter to blastn
    src.fastq.write_fastq2fasta(fq_chunk, query_fpath)
    # Alignm obtain dataframe containing data about alignments, and add 'major' column
    aln_df = _add_major_col(_str2df(src.blast.blast_align(blast_cmd)))

    # If we do remove primers, we need to remove only those
    #   short sequences, which are "major".
    if not primers_lengths is None:
        # Filter short alignments (but only for "major" alignments)
        aln_df = _filter_short_major_alns(
            aln_df,
            min_len_major
        )

        # Select proper function for alignment processing
        #   (this one will check for spurious major alignments,
        #   remove primers and remove short minor amplicons)
        get_aligned_spans = functools.partial(
            _get_aligned_spans_rm_primers,
            primers_lengths=primers_lengths,
            min_len_minor=min_len_minor
        )
    else:
        # Filter short alignments
        aln_df = _filter_short_alns(
            aln_df,
            min_len_major,
            min_len_minor
        )

        # Select proper function for alignment processing
        get_aligned_spans = _get_aligned_spans_keep_primers
    # end if

    result_fastq_records = list()

    for _, fq_record in fq_chunk.items():
        # Select rows containing alignments of current read
        curr_alns = aln_df[aln_df['qseqid'] == fq_record['seq_id']]
        aligned_spans = get_aligned_spans(curr_alns)

        # Write spans
        for aln_span in aligned_spans:

            # Configure variables for output fastq record
            qstart, qend = aln_span[0], aln_span[1]
            # Write 1-based coordinates here
            curr_seq_id = '{}_{}-{}'.format(fq_record['seq_id'], qstart+1, qend)
            curr_seq = fq_record['seq'][qstart : qend]
            curr_qual = fq_record['qual'][qstart : qend]

            result_fastq_records.append({'seq_id':  curr_seq_id,
                                         'seq'   :  curr_seq,
                                         'cmnt'  :  fq_record['cmnt'],
                                         'qual'  :  curr_qual})
        # end for
        with INC_LOCK:
            INC_VAL.value += 1
        # end with

        # Update status bar
        if INC_VAL.value > NEXT_PRINT_NUM.value:
            bar_len = _get_bar_len()
            done_ratio = INC_VAL.value / nreads
            with PRINT_LOCK:
                sys.stdout.write('\r{} - [{}>{}] {}/{} ({}%)'.format(getwt(),
                    '='*int(bar_len*done_ratio),
                    ' '*int(bar_len*(1-done_ratio)), INC_VAL.value, nreads, int(done_ratio*100)))
                sys.stdout.flush()
                NEXT_PRINT_NUM.value += inc_num
            # end with
         # end if
    # end for

    with WRITE_LOCK:
        with open(outfpath, 'a') as outfile:
            for record in result_fastq_records:
                src.fastq.write_fastq_record(record, outfile)
            # end for
        # end with
    # end with

    # Remove temporary query file
    src.filesystem.rm_query_file(query_fpath)
# end def _shredder


def clean_and_shred(
    fq_fpath,
    db_fpath,
    n_thr,
    chunk_size,
    min_len_major,
    min_len_minor,
    primers_lengths):
    # Function organizes analysis of dataframe of alignment data.
    #
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath str;
    # :param db_fpath: path to blast database;
    # :type db_fpath: str;
    # :param n_thr: number of threads;
    # :type n_thr: int;
    # :param chunk_size: size of fastq chunk;
    # :type chunk_size: int;
    # :param min_len_major: minimum length of "major" alignment;
    # :type min_len_major: int;
    # :param min_len_minor: minimum length of "minor" alignment;
    # :type min_len_minor: int;
    # :param primers_lengths: tuple mapping primer's index to it's length.
    #   In this tuple, index of a primer equals index of this primer in source .csv file;
    # :type primers_lengths: tuple<int>;
    #
    # Returns path to output file.

    # Get path to output file.
    outfpath = src.filesystem.make_outfpath(fq_fpath)

    # Empty output file
    with open(outfpath, 'w') as _:
        pass
    # end with

    # Count reads and configure variables for printing status bar
    print('Counting reads...')
    nreads = sum(1 for _ in OPEN_FUNCS[int(fq_fpath.endswith('.gz'))](fq_fpath)) // 4
    bar_len = _get_bar_len()
    inc_num = int(nreads * 0.01)
    print('{} reads.'.format(nreads))

    sys.stdout.write('{} - [{}] 0/{} (0%)'.format(getwt(), ' '*bar_len, nreads))
    sys.stdout.flush()

    # Initialize counters
    global INC_VAL
    global NEXT_PRINT_NUM
    INC_VAL = mp.Value('i', 0)
    NEXT_PRINT_NUM = mp.Value('i', 0)

    # Proceed
    with mp.Pool(n_thr) as pool:
        pool.starmap(_shredder, (
            (
                fq_chunk,
                db_fpath,
                min_len_major,
                min_len_minor,
                nreads,
                inc_num,
                outfpath,
                primers_lengths
            )
            for fq_chunk in src.fastq.fastq_chunks(fq_fpath, chunk_size))
        )
    # end with

    sys.stdout.write('\r{} - [{}] {}/{} (100%)\n'.format(getwt(),
        '='*bar_len, nreads, nreads))
    sys.stdout.flush()

    return outfpath
# end def clean_and_shred
