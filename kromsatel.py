#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.3.b"
# Year, month, day
__last_update_date__ = "2021-01-27"

# |===== Check python interpreter version. =====|

import sys

if sys.version_info.major < 3:
    print( '\nYour python interpreter version is ' + '%d.%d' % (sys.version_info.major,
        sys.version_info.minor) )
    print('   Please, use Python 3.\a')
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        raw_input('Press ENTER to exit:')
    # end if
    sys.exit(1)
else:
    print('\nPython {}.{}.{}\n'.format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
# end if


def platf_depend_exit(exit_code):
    # Function asks to press ENTER press on Windows
    #     and exits after that.

    # :type exit_code: int;

    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def platf_depend_exit


# Firstly check if ve just need to print version or help message
if '-v' in sys.argv[1:] or '--version' in sys.argv[1:] or '-version' in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

if '-h' in sys.argv[1:] or '--help' in sys.argv[1:] or '-help' in sys.argv[1:]:
    print('kromsatel')
    print('Version {}. {} edition.'.format(__version__, __last_update_date__))
    print("Usage:")
    print('  ./kromsatel.py <YOUR_READS> -d <DATABASE_WITH_FRAGMENTS>')
    print('For example:')
    print('  ./kromsatel.py corona_reads.fastq -d fragments-db/nCoV-2019_ref-fragments.fasta')
    platf_depend_exit(0)
# end if

import os
import re
import getopt
import subprocess as sp

import functools
import operator
import io

try:
    import numpy as np
except ImportError:
    print('`numpy` package is not installed')
    print('Please, install it. Example command:')
    print('  pip3 install numpy')
    platf_depend_exit(1)
# end try

try:
    import pandas as pd
except ImportError:
    print('`pandas` package is not installed')
    print('Please, install it. Example command:')
    print('  pip3 install pandas')
    platf_depend_exit(1)
# end try


import gzip
# For opening plain text and gzipped files
OPEN_FUNCS = (open, gzip.open)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode('utf-8').strip()  # format gzipped line
)


from time import time, strftime, gmtime
START_TIME = time() # consider time of importing as start time
def getwt():
    # Function (get work time) returns time HH:MM:SS that has passed from start_time.
    return strftime('%H:%M:%S', gmtime( time() - START_TIME))
# end def getwt


def handle_cl_args():
    # Function handles command line arguments.
    # Returns:
    #  1. List of paths to fastq files to be processed.
    #  2. Path to reference database with fragments.
    #  3. Number of threads to launch.

    # Get arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvd:t:p:', ['help', 'version',
                                                                'db=', 'threads=', 'packet-size=']
        )
    except getopt.GetoptError as opt_err:
        print( str(opt_err) )
        print('See help (`-h` option)')
        platf_depend_exit(2)
    # end try

    # Check and add fastq files:
    fq_fpaths = list()
    is_fastq = lambda f: not re.match(r'.*\.f(ast)?q(\.gz)?$', f) is None

    for arg in args:
        if not os.path.exists(arg):
            print('Error: file `{}` does not exist!'.format(arg))
            platf_depend_exit(1)
        elif not is_fastq(arg):
            print('Error: file `{}` does not look like a fastq file!'.format(arg))
            print('The script detects fastq files by extention and can process gzipped files.')
            print('Permitted extentions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`')
            platf_depend_exit(1)
        else:
            fq_fpaths.append(os.path.abspath(arg))
        # end if
    # end for

    # Check if there are any files to process
    if len(fq_fpaths) == 0:
        print('Usage:')
        print('./kromsatel.py <YOUR_READS> -d <DATABASE_WITH_FRAGMENTS>')
        platf_depend_exit(0)
    # end if

    db_fpath = None
    n_thr = 1
    packet_size = 1000

    # Handle options
    for opt, arg in opts:
        # Get path to database
        if opt in ('-d', '--db'):
            if not os.path.exists( '{}.nhr'.format(arg) ):
                print('Error: database `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                db_fpath = os.path.abspath(arg)
            # end if

        # Handle packet size
        elif opt in ("-p", "--packet-size"):
            try:
                packet_size = int(arg)
                if packet_size <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print("Error: packet_size (-p option) must be positive integer number")
                platf_depend_exit(1)
            # end try

        # Handle number of threads
        elif opt in ('-t', '--threads'):
            try:
                n_thr = int(arg)
                if n_thr < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: number of threads must be positive integer number!')
                print(' And here is your value: `{}`'.format(arg))
                sys.exit(1)
            # end try
            if n_thr > len(os.sched_getaffinity(0)):
                print('''\nWarning! You have specified {} threads to use
      although {} are available.'''.format(n_thr, len(os.sched_getaffinity(0))))
                error = True
                while error:
                    reply = input('''\nPress ENTER to switch to {} threads,
      or enter 'c' to continue with {} threads,
      or enter 'q' to exit:>>'''.format(len(os.sched_getaffinity(0)), n_thr))
                    if reply in ('', 'c', 'q'):
                        error = False
                        if reply == '':
                            n_thr = len(os.sched_getaffinity(0))
                            print('\nNumber of threads switched to {}\n'.format(n_thr))
                        elif reply == 'c':
                            pass
                        elif reply == 'q':
                            sys.exit(0)
                        # end if
                    else:
                        print('\nInvalid reply!\n')
                    # end if
                # end while
            # end if
        # end if
    # end for

    return fq_fpaths, db_fpath, n_thr, packet_size
# end def handle_cl_args


def check_blastn():
    # Check if 'blast+' tookit is installed
    pathdirs = os.environ['PATH'].split(os.pathsep)
    utility = 'blastn'
    utility_found = False
    for directory in pathdirs:
        if os.path.exists(directory) and utility in os.listdir(directory):
            utility_found = True
            break
        # end if
    # end for
    if not utility_found:
        print('  Attention!\n`{}` from BLAST+ toolkit is not installed.'.format(utility))
        print('''If this error still occures although you have installed everything 
-- make sure that this program is added to PATH)''')
        platf_depend_exit(1)
    # end if
# end def check_blastn


def configure_blastn_cmd(query_fpath, db_fpath, n_thr):
    # Function creates command for launching blastn.
    # This command will be always the same, so let's create it once.
    #
    # :param query_fpath: path to query fasta file;
    # :type query_fpath: str;
    # :param db_fpath: path to database;
    # :type db_fpath: str;
    # :param n_thr: number of threads to launch;
    # :type n_thr: int;
    #
    # Returns command which will launch aligning with blastn (str).

    outfmt = '6 qseqid sseqid sstrand length qlen slen qstart qend sstart send'

    # Configure command line
    blast_cmd = "blastn -query {} \
    -db {} \
    -task {} \
    -num_threads {} \
    -evalue 1e-3 \
    -outfmt '{}'".format(query_fpath,
                         db_fpath,
                         'dc-megablast',
                         n_thr,
                         outfmt
    )

    return blast_cmd
# end def configure_blastn_cmd


def disco_align(blast_cmd):
    # Function alignes reads against fragments using Discontiguous Megablast.
    #
    # :param blast_cmd: command to execute;
    # :type blast_cmd: str;
    #
    # Returns tabular string if format "outfmt 6" returned by blastn.

    # Launch blastn
    pipe = sp.Popen(blast_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError while aligning a sequence against amplion database')
        print(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(pipe.returncode)
    # end if

    return stdout_stderr[0].decode('utf-8')
# end def disco_align


def str2df(align_result_str):
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
# end def str2df


def make_outfpath(fq_fpath):
    # Function configures path to output file.
    #
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath: str;
    #
    # Returns path to output file (str).

    try:
        extention = re.search(r'(\.f(ast)?q(\.gz)?)', fq_fpath).group(1)
    except AttributeError as err:
        print( str(err) )
        print('Error -3. Please, contact the developer.')
        platf_depend_exit(-3)
    # end try

    name_with_no_ext = fq_fpath.replace(extention, '')

    return os.path.join(
        os.path.dirname(fq_fpath),
        '{}_cleaned.fastq'.format(name_with_no_ext)
    )
# end def make_outfpath


def form_packet(fastq_file, packet_size, fmt_func):
    """
    Function reads lines from 'fastq_file' and composes a packet of 'packet_size' sequences.

    :param fastq_file: file instance from which to read;
    :type fastq_file: _io.TextIOWrapper or gzip.File;
    :param packet_size: number of sequences to retrive from file;
    :type packet_size: int;
    :param fmt_func: formating function from FORMMATING_FUNCS tuple;
    """

    eof = False
    fq_packet = dict()

    for _ in range(packet_size):

        read_id = fmt_func(fastq_file.readline())

        if read_id == "": # if eof is reached, leave now
            eof = True
            break
        # end if

        fq_packet[read_id] = {
            'seq_id': read_id.partition(' ')[0][1:],
            'seq': fmt_func(fastq_file.readline()),
            'cmnt': fmt_func(fastq_file.readline()),
            'qual': fmt_func(fastq_file.readline())
        }
    # end for

    return fq_packet, eof
# end def form_packet


def fastq_packets(fq_fpath, packet_size):

    how_to_open = OPEN_FUNCS[ fq_fpath.endswith('.gz') ]
    fmt_func = FORMATTING_FUNCS[ fq_fpath.endswith('.gz') ]

    with how_to_open(fq_fpath) as fastq_file:

        # End of file
        eof = False

        while not eof:

            fq_packet, eof = form_packet(fastq_file, packet_size, fmt_func)

            if len(fq_packet) == 0:
                return
            # end if

            yield fq_packet

            if eof:
                return
            # end if
        # end while
    # end with
# end def fastq_packets


def write_fastq2fasta(fq_packet, query_fpath):
    # Function writes fastq-formatted packet to fasta file.
    #
    # :param fq_packet: dictionary of fastq-records;
    # :type fq_packet: dict<str: dict<str: str>>;
    # :param query_fpath: path to query fasta file;
    # :type query_fpath: str;

    with open(query_fpath, 'w') as query_file:
        for fq_record in fq_packet.values():
            query_file.write('>{}\n{}\n'.format(fq_record['seq_id'],fq_record['seq']))
    # end with
# end def


def write_fastq_record(fq_record, outfile):
    # Function writes a fastq record to given file
    #
    # :param fq_record: dictionary of following structure:
    #   {'seq_id': <seq_title>, 'seq': <sequence>, 'cmnt': <comment>, 'qual': <quality_string>};
    # :type fq_record: dict<str: str>;
    # :param outfile: file to which a record will be written;
    # :type outfile: _io.TextIOWrapper;

    outfile.write('@{}\n{}\n{}\n{}\n'.format(fq_record['seq_id'],
            fq_record['seq'], fq_record['cmnt'], fq_record['qual'])
    )
# end def write_fastq_record


# Positions in BLAST+ are 1-based
MAX_EDGE_OFFSET = 10


def set_touch_start(row):
    # Function to be used with pandas.DataFrame.apply function.
    # It adds a column indicating if the alignment in a particular row
    #   "touches" 5'-end of subject sequence.
    #
    # :param row: row of dataframe to which this function if applied;
    # :type row: pandas.Series;

    if row['sstrand']:
        row['touch_start'] = row['sstart'] < 1 + MAX_EDGE_OFFSET
    else:
        row['touch_start'] = row['sstart'] > row['slen'] - MAX_EDGE_OFFSET
    # end if

    return row
# end def get_touch_start


def set_touch_end(row):
    # Function to be used with pandas.DataFrame.apply function.
    # It adds a column indicating if the alignment in a particular row
    #   "touches" 3'-end of subject sequence.
    #
    # :param row: row of dataframe to which this function if applied;
    # :type row: pandas.Series;

    if row['sstrand']:
        row['touch_end'] = row['send'] > row['slen'] - MAX_EDGE_OFFSET
    else:
        row['touch_end'] = row['send'] < 1 + MAX_EDGE_OFFSET
    # end if

    return row
# end def get_touch_start


def set_major(row):
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
# end def get_touch_start


def get_aligned_spans(curr_alns):
    # Function analyses obtained alignments and find "spans" -- subsequences
    #   into which input read should be "shredded".
    # These spans should not overlap
    #
    # :param curr_alns: dataframe containing alignmens data obtained from str2df function;
    # :type curr_alns: pandas.DataFrame;

    # If there are no alignments -- just return empty list.
    if curr_alns.empty:
        return []
    # end if

    # Add columns indicating if alignments "touch" 5'- and 3'-ends of subjects
    # And a column indicating if subject sequence is a major fragment
    curr_alns.insert(loc=curr_alns.shape[1], column='touch_start', value=np.repeat(None, curr_alns.shape[0]))
    curr_alns.insert(loc=curr_alns.shape[1], column='touch_end', value=np.repeat(None, curr_alns.shape[0]))
    curr_alns.insert(loc=curr_alns.shape[1], column='major', value=np.repeat(False, curr_alns.shape[0]))

    # Fill these columns
    curr_alns = curr_alns.apply(set_touch_start, axis=1)
    curr_alns = curr_alns.apply(set_touch_end, axis=1)
    curr_alns = curr_alns.apply(set_major, axis=1)

    # Extract "touching" alignments
    full_span_alns = curr_alns[(curr_alns['touch_start']) & (curr_alns['touch_end'])] # from 5' to 3'
    one_side_alns = curr_alns[operator.xor(curr_alns['touch_start'], curr_alns['touch_end'])] # from 5' of from 3' (with a break somewhere in between)

    chain_of_alns = (
        full_span_alns[full_span_alns['major'] == True], # firstly check major full alignments
        one_side_alns[one_side_alns['major'] == True], # then major but not full
        full_span_alns[full_span_alns['major'] == False], # then minor full
        one_side_alns[one_side_alns['major'] == False], # and finally minor but not full
    )

    del full_span_alns
    del one_side_alns

    # List of result spans
    aligned_spans = list()

    do_overlaps = lambda span: span[0] <= aln['qstart']-1 <= span[1]\
                            or span[0] <= aln['qend']     <= span[1]

    for aln_collection in chain_of_alns:
        for i in range(aln_collection.shape[0]):
            aln = aln_collection.iloc[i, [6, 7]] # get current alignment

            overlap_result = functools.reduce(operator.or_,
                map(do_overlaps, aligned_spans),
                False
            )

            if not overlap_result:
                aligned_spans.append( (aln['qstart']-1, aln['qend']) )
            # end if
        # end for
    # end for

    return aligned_spans
# end def get_aligned_spans


def rm_query_file(query_fpath):
    # Function removes temporary query file.
    #
    # :param query_fpath: path to query file;
    # :type query_fpath: str;

    try:
        os.unlink(query_fpath)
    except OSError as oserr:
        print('Warning: Cannot remove temporary file `{}`'.format(query_fpath))
        print( str(oserr) )
    # end try
# end rm_query_file


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


def clean_and_shred(fq_fpath, db_fpath, n_thr, packet_size):
    # Function organizes analysis of dataframe of alignment data.
    #
    # :param aln_df: dataframe containing alignment data to be analysed;
    # :type aln_df: pandas.DatFrame;
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath str;
    #
    # Returns path to output file.

    # Get path to output file.
    outfpath = make_outfpath(fq_fpath)
    if sys.platform.startswith('win'):
        query_fpath = 'kromsatel_query_{}.fasta'.format(os.getpid())
    else:
        query_fpath = '/tmp/kromsatel_query_{}.fasta'.format(os.getpid())
    # end if

    blast_cmd = configure_blastn_cmd(query_fpath, db_fpath, n_thr)

    # Count reads and configure variables for printing status bar
    nreads = sum(1 for _ in OPEN_FUNCS[int(fq_fpath.endswith('.gz'))](fq_fpath)) // 4
    bar_len = _get_bar_len()
    next_print_num = int(nreads * 0.01)
    inc_num = next_print_num
    i = 0

    sys.stdout.write('{} - [{}] 0/{} (0%)'.format(getwt(), ' '*bar_len, nreads))
    sys.stdout.flush()

    with open(outfpath, 'w') as outfile:

        for fq_packet in fastq_packets(fq_fpath, packet_size):

            # Convert input fastq file to fasta format in order to pass the latter to blastn
            write_fastq2fasta(fq_packet, query_fpath)
            # Align and obtain dataframe containing data about alignments
            aln_df = str2df(disco_align(blast_cmd))

            for _, fq_record in fq_packet.items():
                # Select rows containing alignments of current read
                curr_alns = aln_df[aln_df['qseqid'] == fq_record['seq_id']]
                aligned_spans = get_aligned_spans(curr_alns)

                # Write spans
                for aln_span in aligned_spans:

                    # Configure variables for output fastq record
                    qstart, qend = aln_span[0], aln_span[1]
                    curr_seq_id = '{}_{}-{}'.format(fq_record['seq_id'], qstart+1, qend) # write 1-based coordinates here
                    curr_seq = fq_record['seq'][qstart : qend]
                    curr_qual = fq_record['qual'][qstart : qend]

                    write_fastq_record({'seq_id':  curr_seq_id,
                                        'seq'   :  curr_seq,
                                        'cmnt'  :  fq_record['cmnt'],
                                        'qual'  :  curr_qual},
                                       outfile)
                # end for
                i += 1

                # Update status bar
                if i > next_print_num:
                    bar_len = _get_bar_len()
                    done_ratio = i / nreads
                    sys.stdout.write('\r{} - [{}>{}] {}/{} ({}%)'.format(getwt(),
                        '='*int(bar_len*done_ratio),
                        ' '*int(bar_len*(1-done_ratio)), i, nreads, int(done_ratio*100)))
                    sys.stdout.flush()
                    next_print_num += inc_num
                 # end if
            # end for
        # end for
    # end with

    sys.stdout.write('\r{} - [{}] {}/{} (100%)\n'.format(getwt(),
        '='*bar_len, nreads, nreads))
    sys.stdout.flush()

    # Remove temporary query file
    rm_query_file(query_fpath)

    return outfpath
# end def clean_and_shred


def main():
    # Handle command line arguments
    fq_fpaths, db_fpath, n_thr, packet_size = handle_cl_args()

    # Check if blastn is installed
    check_blastn()

    print('{} - Start.'.format(getwt()))

    for fq_fpath in fq_fpaths:
        print('{} - Processing file `{}`'.format(getwt(), fq_fpath))

        # Run cleaning
        outfpath = clean_and_shred(fq_fpath, db_fpath, n_thr, packet_size)

        print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
        print('Output file: `{}`'.format(outfpath))
        print('-------')
    # end for

    print('{} - Completed.'.format(getwt()))
# end def main


if __name__ == '__main__':
    main()
# end if
