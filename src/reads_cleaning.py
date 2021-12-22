
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



class ReadsCleaner:

    def __init__(self, args):
        self.args = args

        num_reads_total = self._count_reads()

        self.primer_scheme = prm.PrimerScheme(self.args)
        self.progress = Progress(num_reads_total)
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
        print()

        forward_chunk = reads_chunk[0]
        forward_alignments = parse_alignments(
            src.blast.blast_align(forward_chunk, self.args)
        )

        reverse_chunk = reads_chunk[1]
        reverse_alignments = parse_alignments(
            src.blast.blast_align(reverse_chunk, self.args)
        )

        # print()
        # print(len(forward_alignments))
        # print(len(reverse_alignments))
        # f_ids = map(lambda x: x.partition('__<SPACE>__')[0], forward_alignments.keys())
        # r_ids = map(lambda x: x.partition('__<SPACE>__')[0], reverse_alignments.keys())
        # print(
        #     set(f_ids) == set(r_ids)
        # )

        # for fn, rn in zip(forward_alignments.keys(), reverse_alignments.keys()):
        #     fn_p = fn.partition('__<SPACE>__')[0]
        #     rn_p = rn.partition('__<SPACE>__')[0]
        #     if fn_p != rn_p:
        #         print(fn)



        for forward_read, reverse_read in zip(*reads_chunk):
            forward_alignment = forward_alignments[forward_read['seq_id']]
            reverse_alignment = reverse_alignments[reverse_read['seq_id']]

            if not (forward_alignment is None or reverse_alignment is None):

                if forward_alignment.align_strand_plus:
                    forward_start_coord = forward_alignment.ref_start
                    forward_end_coord   = forward_alignment.ref_end
                else:
                    forward_start_coord = forward_alignment.ref_end
                    forward_end_coord   = forward_alignment.ref_start
                # end if
                if reverse_alignment.align_strand_plus:
                    reverse_start_coord = reverse_alignment.ref_start
                    reverse_end_coord   = reverse_alignment.ref_end
                else:
                    reverse_start_coord = reverse_alignment.ref_end
                    reverse_end_coord   = reverse_alignment.ref_start
                # end if


                MAJOR, MINOR, ABNORMAL = range(3)
                classification = ABNORMAL
                FIXED_SHRINK_LEN = 27
                forward_end_needs_shrink, reverse_end_needs_shrink = True, True

                forward_start_primer_num, reverse_start_primer_num = None, None
                forward_end_primer_num, reverse_end_primer_num = None, None

                forward_left = forward_alignment.align_strand_plus
                reverse_left = reverse_alignment.align_strand_plus

                forward_start_primer_num, _ = \
                    self.primer_scheme.find_primer_by_coord(
                        forward_start_coord
                    )

                forward_start_primer_found = not forward_start_primer_num is None
                reverse_left = not forward_left

                print(forward_alignment)
                print(reverse_alignment)

                if forward_start_primer_found:
                    reverse_start_major = \
                        self.primer_scheme.check_coord_within_primer(
                            reverse_start_coord,
                            forward_start_primer_num,
                            left=reverse_left
                        )
                    if reverse_start_major:
                        classification = MAJOR
                        reverse_start_primer_num = forward_start_primer_num
                    else:
                        minor_pair_primer_num = forward_start_primer_num + (-1 if forward_left else 1)
                        reverse_start_minor = \
                            self.primer_scheme.check_coord_within_primer(
                                reverse_start_coord,
                                minor_pair_primer_num,
                                left=reverse_left
                            )
                        if reverse_start_minor:
                            classification = MINOR
                            reverse_start_primer_num = minor_pair_primer_num
                        # end if
                    # end if
                else:
                    reverse_start_primer_num, _ = \
                        self.primer_scheme.find_primer_by_coord(
                            reverse_start_coord
                        )
                # end if

                if not forward_start_primer_num is None:
                    reverse_end_needs_shrink = self.primer_scheme.check_coord_within_primer(
                        reverse_end_coord,
                        forward_start_primer_num,
                        left=forward_left
                    )
                    if reverse_end_needs_shrink:
                        reverse_end_primer_num = forward_start_primer_num
                    # end if
                # end if
                if not reverse_start_primer_num is None:
                    forward_end_needs_shrink = self.primer_scheme.check_coord_within_primer(
                        forward_end_coord,
                        reverse_start_primer_num,
                        left=reverse_left
                    )
                    if forward_end_needs_shrink:
                        forward_end_primer_num = reverse_start_primer_num
                    # end if
                # end if


                print('Classification: {}'.format(classification))
                print('Forward start primer', forward_start_primer_num, forward_left)
                print('Reverse start primer', reverse_start_primer_num, reverse_left)
                print('Forward end needs shrink: {} ({})'.format(forward_end_needs_shrink, forward_end_primer_num))
                print('Reverse end needs shrink: {} ({})'.format(reverse_end_needs_shrink, reverse_end_primer_num))

                if classification == 0:
                    reverse_start_primer_num = forward_start_primer_num
                elif classification == 1:
                    reverse_start_primer_num = minor_pair_primer_num
                # end if

                if not forward_start_primer_num is None:
                    forward_alignment = self._shrink_start(forward_alignment, forward_start_primer_num, forward_left)
                # end if
                if not reverse_start_primer_num is None:
                    reverse_alignment = self._shrink_start(reverse_alignment, reverse_start_primer_num, reverse_left)
                # end if
                if not forward_end_primer_num is None:
                    forward_alignment = self._shrink_end(forward_alignment, forward_end_primer_num, (not forward_left))
                # end if
                if not reverse_end_primer_num is None:
                    reverse_alignment = self._shrink_end(reverse_alignment, reverse_end_primer_num, (not reverse_left))
                # end if
                print('After shrinking:')
                print(forward_alignment)
                print(reverse_alignment)
                print()


                # forward_start_primer = self.primer_scheme.find_primer_by_coord(forward_alignment.ref_start)
                # forward_end_primer = self.primer_scheme.find_primer_by_coord(forward_alignment.ref_end)

                # reverse_start_primer = self.primer_scheme.find_primer_by_coord(reverse_alignment.ref_start)
                # reverse_end_primer = self.primer_scheme.find_primer_by_coord(reverse_alignment.ref_end)

                # print(forward_alignment)
                # print(forward_start_primer, forward_end_primer)
                # print(reverse_alignment)
                # print(reverse_start_primer, reverse_end_primer)
                # print()
            # end if
        # end for
    # end def


    def _shrink_start(self, alignment, primer_num, left=True):
        if left:
            primer = self.primer_scheme.primer_pairs[primer_num].left_primer
        else:
            primer = self.primer_scheme.primer_pairs[primer_num].right_primer
        # end if

        if alignment.align_strand_plus:
            primer_len_in_read = primer.end - alignment.ref_start + 1
            alignment.ref_start += primer_len_in_read
            alignment.query_start += primer_len_in_read
        else:
            primer_len_in_read = alignment.ref_end - primer.start + 1
            alignment.ref_end -= primer_len_in_read
            alignment.query_end -= primer_len_in_read
        # end if

        return alignment
    # end def

    def _shrink_end(self, alignment, primer_num, left=True):
        if left:
            primer = self.primer_scheme.primer_pairs[primer_num].left_primer
        else:
            primer = self.primer_scheme.primer_pairs[primer_num].right_primer
        # end if

        if alignment.align_strand_plus:
            primer_len_in_read = alignment.ref_end - primer.start + 1
            alignment.ref_end -= primer_len_in_read
            alignment.query_end -= primer_len_in_read
        else:
            primer_len_in_read = primer.end - alignment.ref_start + 1
            alignment.ref_start += primer_len_in_read
            alignment.query_start += primer_len_in_read
        # end if

        return alignment
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

# end class



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
        ratio_done = curr_num_done_reads / self.NUM_READS_TOTAL
        progress_line_len = round(bar_len * ratio_done)


        global print_lock
        with print_lock:
            sys.stdout.write(
                '\r{} - [{}{}] {}/{} ({}%)'.format(
                    getwt(),
                    '=' * progress_line_len,
                    ' ' * (bar_len - progress_line_len),
                    curr_num_done_reads,
                    self.NUM_READS_TOTAL,
                    round(ratio_done * 100)
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