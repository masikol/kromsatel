
import os
import sys
import functools
import multiprocessing as mp

import src.fastq
import src.primers as prm
from src.printing import getwt
from src.alignment import parse_alignments, Alignment


# This damned stuff should be global, because
#   "Synchronized objects should only be shared between processes through inheritance"
num_done_reads = mp.Value('i', 0)
next_report_num = mp.Value('i', 0)
output_lock = mp.Lock()
print_lock = mp.Lock()
status_update_lock = mp.Lock()


class Progress:

    def __init__(self, num_reads_total):
        self.NUM_READS_TOTAL = num_reads_total
        self._REPORT_DELAY = round(self.NUM_READS_TOTAL * 0.01)
        self._DEFAULT_STATUS_BAR_LEN = 40

        global next_report_num
        next_report_num.value = self._REPORT_DELAY
    # end def __init__

    def get_num_done_reads(self):
        global num_done_reads
        return num_done_reads.value
    # end def get_num_done_reads

    def get_next_report_num(self):
        global next_report_num
        return next_report_num.value
    # end def get_next_report_num

    def increment_done(self, increment=1):
        global num_done_reads
        num_done_reads.value = num_done_reads.value + increment
    # end def increment_done

    def increment_next_report(self):
        global next_report_num
        next_report_num = next_report_num.value + self._REPORT_DELAY
    # end def increment_next_report

    def print_status_bar(self):

        curr_num_done_reads = self.get_num_done_reads()

        bar_len = self._get_status_bar_len()
        percent_done = round(curr_num_done_reads / self.NUM_READS_TOTAL * 100)

        global print_lock
        with print_lock:
            sys.stdout.write(
                '\r{} - [{}] {}/{} ({}%)'.format(
                    getwt(),
                    '=' * bar_len,
                    curr_num_done_reads,
                    self.NUM_READS_TOTAL,
                    percent_done
                )
            )
            sys.stdout.flush()
        # end with
    # end def print_status_bar

    def _get_status_bar_len(self):
        try:
            bar_len = int(os.get_terminal_size().columns * 0.40)
        except OSError:
            bar_len = self._DEFAULT_STATUS_BAR_LEN
        # end try
        return bar_len
    # end _get_status_bar_len
# end class Progress


class ReadsCleaner:

    def __init__(self, args):
        self.args = args

        num_reads_total = self._count_reads()

        self.progress = Progress(num_reads_total)
        self.primer_scheme = prm.PrimerScheme(self.args)

        self._get_chunk_interator = self._choose_chunk_interator()
    # end def __init__


    def clean_reads(self):

        reads_chunks = self._choose_fastq_chunks_func()

        self.progress.print_status_bar()

        # Proceed
        with mp.Pool(self.args['n_thr']) as pool:
            pool.starmap(
                self._clean_reads_chunk_paired,
                (
                    (reads_chunk,) for reads_chunk in reads_chunks()
                )
            )
        # end with

        pool.close()
        pool.join()

        self.progress.print_status_bar()
        print()
    # end def clean_reads


    def _clean_reads_chunk_paired(self, reads_chunk):

        chunk_interator = self._get_chunk_interator(reads_chunk)

        forward_chunk = reads_chunk[0]
        forward_alignments = parse_alignments(
            src.blast.blast_align(forward_chunk, self.args)
        )

        reverse_chunk = reads_chunk[1]
        reverse_alignments = parse_alignments(
            src.blast.blast_align(reverse_chunk, self.args)
        )

        print()
        print(len(forward_alignments))
        print(len(reverse_alignments))
        f_ids = map(lambda x: x.partition('__<SPACE>__')[0], forward_alignments.keys())
        r_ids = map(lambda x: x.partition('__<SPACE>__')[0], reverse_alignments.keys())
        print(
            set(f_ids) == set(r_ids)
        )

        for fn, rn in zip(forward_alignments.keys(), reverse_alignments.keys()):
            fn_p = fn.partition('__<SPACE>__')[0]
            rn_p = rn.partition('__<SPACE>__')[0]
            if fn_p != rn_p:
                print(fn)
                


        # for read_or_pair in chunk_interator:
        #     self._classify_read_pair(read_or_pair)
        # # end for
    # end def

    def _count_reads(self):
        print('{} - Counting reads...'.format(getwt()))
        if self.args['paired_mode']:
            num_reads_total = src.fastq.count_reads(self.args['reads_R1'])
        else:
            num_reads_total = src.fastq.count_reads(self.args['reads_unpaired'])
        # end if
        print('{} - {} reads.'.format(getwt(), num_reads_total))
        return num_reads_total
    # end def _count_reads


    def _choose_fastq_chunks_func(self):
        if self.args['paired_mode']:
            fastq_chunks = functools.partial(
                src.fastq.fastq_chunks_paired,
                forward_read_fpath=self.args['reads_R1'],
                reverse_read_fpath=self.args['reads_R2'],
                chunk_size=self.args['chunk_size']
            )
        else:
            fastq_chunks = functools.partial(
                src.fastq.fastq_chunks_unpaired,
                fq_fpath=self.args['reads_unpaired'],
                chunk_size=self.args['chunk_size']
            )
        # end if
        return fastq_chunks
    # end def _choose_fastq_chunks_func


    def _choose_chunk_interator(self):

        if self.args['paired_mode']:
            get_chunk_interator = self._get_paired_chunk_iterator
        else:
            get_chunk_interator = self._get_unpaired_chunk_iterator
        # end if
        return get_chunk_interator
    # end def _choose_chunk_interator


    def _get_paired_chunk_iterator(self, reads_chunk):
        return zip(*reads_chunk)
    # end def


    def _get_unpaired_chunk_iterator(self, reads_chunk):
        return reads_chunk
    # end def


    def _classify_read_pair(self, read_pair):

        forward_read = read_pair[0]
        reverse_read = read_pair[1]

        forward_left = True

        forward_primer_class = self.primer_scheme.find_left_primer(
            forward_read
        )

        if forward_primer_class is None:
            forward_left = False
        # end if

        if forward_left:
            reverse_primer_class = self.primer_scheme.find_right_primer(
                reverse_read
            )
        else:
            forward_primer_class = self.primer_scheme.find_right_primer(
                forward_read
            )
            reverse_primer_class = self.primer_scheme.find_left_primer(
                reverse_read
            )
        # end if

        print(forward_read['seq_id'])
        print(forward_primer_class, reverse_primer_class, forward_left)
    # end def


    # def _classify_read(self, read):


    # # end def

# end class
