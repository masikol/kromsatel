
import os

import src.filesystem as fs


SPACE_HOLDER = '__<SPACE>__'


class FastqRecord:

    def __init__(self, header, seq, comment, quality_str):
        self.header = header
        self.seq = seq
        self.comment = comment
        self.quality_str = quality_str
    # end def

    def get_seqid(self):
        return self.header.partition(SPACE_HOLDER)[0]
    # end def

    def get_copy(self):
        return FastqRecord(
            self.header, self.seq,
            self.comment, self.quality_str
        )
    # end def

    def __len__(self):
        return len(self.seq)
    # end def
# end class


def count_reads(file_path):
    is_gzipped = file_path.endswith('.gz')
    return sum(
        1 for _ in fs.open_file_may_by_gzipped(file_path)
    ) // 4
# end def


def form_chunk(fastq_file, chunk_size):

    eof = False
    fq_chunk = [None] * chunk_size

    for i in range(chunk_size):

        header = fastq_file.readline().strip()

        if header == '': # if eof is reached, terminate reading
            eof = True
            break
        # end if

        formatted_header = header[1:].replace(' ', SPACE_HOLDER)
        seq         = fastq_file.readline().strip().upper()
        comment     = fastq_file.readline().strip()
        quality_str = fastq_file.readline().strip()

        fq_chunk[i] = FastqRecord(
            formatted_header,
            seq,
            comment,
            quality_str
        )
    # end for

    not_none = lambda x: not x is None

    return tuple(filter(not_none, fq_chunk)), eof
# end def


def fastq_chunks_unpaired(fq_fpath, chunk_size):

    with fs.open_file_may_by_gzipped(fq_fpath, 'rt') as fastq_file:

        eof = False # end of file

        while not eof:

            fq_chunk, eof = form_chunk(fastq_file, chunk_size)

            if len(fq_chunk) == 0:
                return
            # end if

            yield fq_chunk

            if eof:
                return
            # end if
        # end while
    # end with
# end def


def fastq_chunks_paired(forward_read_fpath, reverse_read_fpath, chunk_size):

    with fs.open_file_may_by_gzipped(forward_read_fpath) as forward_file, \
         fs.open_file_may_by_gzipped(reverse_read_fpath) as reverse_file:

        eof = False

        while not eof:

            forward_chunk, f_eof = form_chunk(forward_file, chunk_size)
            reverse_chunk, r_eof = form_chunk(reverse_file, chunk_size)

            if len(forward_chunk) == 0 or len(reverse_chunk) == 0:
                return
            # end if

            yield (forward_chunk, reverse_chunk)

            eof = f_eof or r_eof
        # end while
    # end with
# end def


def write_fastq2fasta(reads_chunk, query_fpath):

    with open(query_fpath, 'wt') as query_file:
        for fq_record in reads_chunk:
            query_file.write('>{}\n{}\n'.format(fq_record.header, fq_record.seq))
    # end with
# end def


def write_fastq_record(fq_record, outfile):
    # `outfile` should be opened for appending
    outfile.write('@{}\n{}\n{}\n{}\n'.format(
        fq_record.header.replace(SPACE_HOLDER, ' '),
        fq_record.seq, fq_record.comment, fq_record.quality_str
        )
    )
# end def
