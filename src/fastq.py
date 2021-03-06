# -*- coding: utf-8 -*-

from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS


def form_chunk(fastq_file, chunk_size, fmt_func):
    # Function reads lines from 'fastq_file' and composes a chunk of 'chunk_size' sequences.
    #
    # :param fastq_file: file instance from which to read;
    # :type fastq_file: _io.TextIOWrapper or gzip.File;
    # :param chunk_size: number of sequences to retrive from file;
    # :type chunk_size: int;
    # :param fmt_func: formating function from FORMMATING_FUNCS tuple;

    eof = False
    fq_chunk = dict()

    for _ in range(chunk_size):

        read_id = fmt_func(fastq_file.readline())

        if read_id == "": # if eof is reached, leave now
            eof = True
            break
        # end if

        fq_chunk[read_id] = {
            'seq_id': read_id.partition(' ')[0][1:],
            'seq': fmt_func(fastq_file.readline()),
            'cmnt': fmt_func(fastq_file.readline()),
            'qual': fmt_func(fastq_file.readline())
        }
    # end for

    return fq_chunk, eof
# end def form_chunk


def fastq_chunks(fq_fpath, chunk_size):

    how_to_open = OPEN_FUNCS[ fq_fpath.endswith('.gz') ]
    fmt_func = FORMATTING_FUNCS[ fq_fpath.endswith('.gz') ]

    with how_to_open(fq_fpath) as fastq_file:

        # End of file
        eof = False

        while not eof:

            fq_chunk, eof = form_chunk(fastq_file, chunk_size, fmt_func)

            if len(fq_chunk) == 0:
                return
            # end if

            yield fq_chunk

            if eof:
                return
            # end if
        # end while
    # end with
# end def fastq_chunks


def write_fastq2fasta(fq_chunk, query_fpath):
    # Function writes fastq-formatted chunk to fasta file.
    #
    # :param fq_chunk: dictionary of fastq-records;
    # :type fq_chunk: dict<str: dict<str: str>>;
    # :param query_fpath: path to query fasta file;
    # :type query_fpath: str;

    with open(query_fpath, 'w') as query_file:
        for fq_record in fq_chunk.values():
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
