
import os

from src.filesystem import OPEN_FUNCS, FORMATTING_FUNCS


SPACE_HOLDER = '__<SPACE>__'


def count_reads(file_path):
    is_gzipped = file_path.endswith('.gz')
    return sum(1 for _ in OPEN_FUNCS[int(is_gzipped)](file_path)) // 4
# end def count_reads


def form_chunk(fastq_file, chunk_size, fmt_func):
    # Function reads lines from 'fastq_file' and composes a chunk of 'chunk_size' sequences.
    #
    # :param fastq_file: file instance from which to read;
    # :type fastq_file: _io.TextIOWrapper or gzip.File;
    # :param chunk_size: number of sequences to retrive from file;
    # :type chunk_size: int;
    # :param fmt_func: formating function from FORMMATING_FUNCS tuple;

    eof = False
    fq_chunk = [None] * chunk_size

    for i in range(chunk_size):

        read_id = fmt_func(fastq_file.readline())

        if read_id == "": # if eof is reached, leave now
            eof = True
            break
        # end if

        fq_chunk[i] = {
            'seq_id': read_id[1:].replace(' ', SPACE_HOLDER),
            'seq': fmt_func(fastq_file.readline()),
            'cmnt': fmt_func(fastq_file.readline()),
            'qual': fmt_func(fastq_file.readline())
        }
    # end for

    not_none = lambda x: not x is None

    print('Chunk formed', os.getpid())

    return tuple(filter(not_none, fq_chunk)), eof
# end def form_chunk


def fastq_chunks_unpaired(fq_fpath, chunk_size):

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
# end def fastq_chunks_unpaired


def fastq_chunks_paired(forward_read_fpath, reverse_read_fpath, chunk_size):

    how_to_open_forward = OPEN_FUNCS[ forward_read_fpath.endswith('.gz') ]
    fmt_func_forward = FORMATTING_FUNCS[ forward_read_fpath.endswith('.gz') ]
    how_to_open_reverse = OPEN_FUNCS[ reverse_read_fpath.endswith('.gz') ]
    fmt_func_reverse = FORMATTING_FUNCS[ reverse_read_fpath.endswith('.gz') ]

    with how_to_open_forward(forward_read_fpath) as forward_file, \
         how_to_open_reverse(reverse_read_fpath) as reverse_file:

        # End of file
        eof = False

        while not eof:

            forward_chunk, f_eof = form_chunk(forward_file, chunk_size, fmt_func_forward)
            reverse_chunk, r_eof = form_chunk(reverse_file, chunk_size, fmt_func_reverse)

            if len(forward_chunk) == 0 or len(reverse_chunk) == 0:
                return
            # end if

            yield (forward_chunk, reverse_chunk)

            if f_eof or r_eof:
                return
            # end if
        # end while
    # end with
# end def fastq_chunks_paired


def write_fastq2fasta(reads_chunk, query_fpath):

    with open(query_fpath, 'w') as query_file:
        # for fq_record in fq_chunk.values():
        for fq_record in reads_chunk:
            query_file.write('>{}\n{}\n'.format(fq_record['seq_id'],fq_record['seq']))
    # end with
# end def


def write_fastq_record(fq_record, outfile):
    # `outfile` should be opened for appending
    outfile.write('@{}\n{}\n{}\n{}\n'.format(
        fq_record['seq_id'].replace(SPACE_HOLDER, ' '),
        fq_record['seq'], fq_record['cmnt'], fq_record['qual']
        )
    )
# end def
