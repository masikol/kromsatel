
import os
import sys
import functools
import multiprocessing as mp

import src.fastq
import src.primers as prm
from src.printing import getwt
from src.alignment import parse_alignments, Alignment

from src.binning import PairedBinner
from src.binning import MAJOR, MINOR, ABNORMAL


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

        self.MIN_LEN = args['min_len']
        self.threads = args['n_thr']
        self.FIXED_SHRINK_LEN = 27 # bp

        num_reads_total = self._count_reads()

        self.primer_scheme = prm.PrimerScheme(args)
        self.progress = Progress(num_reads_total)
    # end def __init__


    def clean_reads(self):

        reads_chunks = self._choose_fastq_chunks_func()

        self.progress.print_status_bar()

        # Proceed
        with mp.Pool(self.threads) as pool:
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
        # print()
    # end def clean_reads


    def _clean_reads_chunk_paired(self, reads_chunk):
        # print()

        forward_chunk = reads_chunk[0]
        forward_alignments = parse_alignments(
            src.blast.blast_align(forward_chunk, self.args)
        )

        reverse_chunk = reads_chunk[1]
        reverse_alignments = parse_alignments(
            src.blast.blast_align(reverse_chunk, self.args)
        )

        binner = PairedBinner(self.args['output'])

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

                classification = ABNORMAL
                crop_forward_end, crop_reverse_end = True, True

                forward_start_primer_num, reverse_start_primer_num = None, None
                forward_end_primer_num, reverse_end_primer_num = None, None

                forward_left = forward_alignment.align_strand_plus
                reverse_left = reverse_alignment.align_strand_plus

                if forward_left:
                    forward_start_primer_num = \
                        self.primer_scheme.find_left_primer_by_coord(
                            forward_start_coord
                        )
                else:
                    forward_start_primer_num = \
                        self.primer_scheme.find_right_primer_by_coord(
                            forward_start_coord
                        )
                # end if

                forward_start_primer_found = not forward_start_primer_num is None
                reverse_left = not forward_left

                # print(forward_alignment)
                # print(reverse_alignment)

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
                    if reverse_left:
                        reverse_start_primer_num = \
                            self.primer_scheme.find_left_primer_by_coord(
                                reverse_start_coord
                            )
                    else:
                        reverse_start_primer_num = \
                            self.primer_scheme.find_right_primer_by_coord(
                                reverse_start_coord
                            )
                    # end if
                # end if

                if not forward_start_primer_num is None:
                    crop_reverse_end = self.primer_scheme.check_coord_within_primer(
                        reverse_end_coord,
                        forward_start_primer_num,
                        left=forward_left
                    )
                    if crop_reverse_end:
                        reverse_end_primer_num = forward_start_primer_num
                    # end if
                # end if
                if not reverse_start_primer_num is None:
                    crop_forward_end = self.primer_scheme.check_coord_within_primer(
                        forward_end_coord,
                        reverse_start_primer_num,
                        left=reverse_left
                    )
                    if crop_forward_end:
                        forward_end_primer_num = reverse_start_primer_num
                    # end if
                # end if

                # print('Classification: {}'.format(classification))
                # print('Forward start primer', forward_start_primer_num, forward_left)
                # print('Reverse start primer', reverse_start_primer_num, reverse_left)
                # print('Forward end needs crop: {} ({})'.format(crop_forward_end, forward_end_primer_num))
                # print('Reverse end needs crop: {} ({})'.format(crop_reverse_end, reverse_end_primer_num))

                if classification == MAJOR:
                    reverse_start_primer_num = forward_start_primer_num
                elif classification == MINOR:
                    reverse_start_primer_num = minor_pair_primer_num
                # end if

                if not forward_start_primer_num is None:
                    forward_alignment = self._trim_start_primer(forward_alignment, forward_start_primer_num, forward_left)
                else:
                    forward_alignment = self._crop_start(forward_alignment)
                # end if
                if not forward_end_primer_num is None:
                    forward_alignment = self._trim_end_primer(forward_alignment, forward_end_primer_num, (not forward_left))
                else:
                    if crop_forward_end:
                        forward_alignment = self._crop_end(forward_alignment)
                    # end if
                # end if

                if not reverse_start_primer_num is None:
                    reverse_alignment = self._trim_start_primer(reverse_alignment, reverse_start_primer_num, reverse_left)
                else:
                    reverse_alignment = self._crop_start(reverse_alignment)
                # end if
                if not reverse_end_primer_num is None:
                    reverse_alignment = self._trim_end_primer(reverse_alignment, reverse_end_primer_num, (not reverse_left))
                else:
                    if crop_reverse_end:
                        reverse_alignment = self._crop_end(reverse_alignment)
                    # end if
                # end if

                # print('After trimming and cropping:')
                # print(forward_alignment)
                # print(reverse_alignment)

                trimmed_forward_align_len = forward_alignment.get_align_len()
                trimmed_reverse_align_len = reverse_alignment.get_align_len()
                forward_survives = trimmed_forward_align_len >= self.MIN_LEN
                reverse_survives = trimmed_reverse_align_len >= self.MIN_LEN
                # print('Alignment lengths: F:{}, R:{}'.format(trimmed_forward_align_len, trimmed_reverse_align_len))
                # print('Survivors: F:{}, R:{}'.format(forward_survives, reverse_survives))

                if forward_survives:
                    forward_read = self._trim_read(forward_read, forward_alignment)
                # end if
                if reverse_survives:
                    reverse_read = self._trim_read(reverse_read, reverse_alignment)
                # end if

                # Binning
                if forward_survives and reverse_survives:
                    if classification == MAJOR:
                        binner.add_major_pair(forward_read, reverse_read)
                    elif classification == MINOR:
                        binner.add_minor_pair(forward_read, reverse_read)
                    else:
                        binner.add_abnormal_pair(forward_read, reverse_read)
                    # end if
                elif forward_survives and not reverse_survives:
                    binner.add_forward_unpaired_read(forward_read)
                elif not forward_survives and reverse_survives:
                    binner.add_reverse_unpaired_read(reverse_read)
                # end if
            # end if
            # print('\n')
        # end for

        # print()
        # print('Major:')
        # print(binner.major_forward_reads)
        # print(binner.major_reverse_reads)
        # print('Minor:')
        # print(binner.minor_forward_reads)
        # print(binner.minor_reverse_reads)
        # print('Abnormal:')
        # print(binner.abnormal_forward_reads)
        # print(binner.abnormal_reverse_reads)
        # print('Unpaired:')
        # print(binner.unpaired_forward_reads)
        # print(binner.unpaired_reverse_reads)

        with output_lock:
            binner.write_binned_reads()
        # end with
        with status_update_lock:
            prev_next_value = self.progress.get_next_report_num()
            self.progress.increment_done(len(forward_chunk))
            self.progress.increment_next_report()
        # end with
        if self.progress.get_next_report_num() != prev_next_value:
            with print_lock:
                self.progress.print_status_bar()
            # end with
        # end if
    # end def

    def _trim_read(self, read, alignment):
        new_start, new_end = alignment.query_start, alignment.query_end+1
        read['seq']  = read['seq'] [new_start : new_end]
        read['qual'] = read['qual'][new_start : new_end]
        return read
    # end def


    def _crop_start(self, alignment):
        if alignment.align_strand_plus:
            alignment.ref_start += self.FIXED_SHRINK_LEN
            alignment.query_start += self.FIXED_SHRINK_LEN
        else:
            alignment.ref_end -= self.FIXED_SHRINK_LEN
            alignment.query_start += self.FIXED_SHRINK_LEN
        # end if
        return alignment
    # end def

    def _crop_end(self, alignment):
        if alignment.align_strand_plus:
            alignment.ref_end -= self.FIXED_SHRINK_LEN
            alignment.query_end -= self.FIXED_SHRINK_LEN
        else:
            alignment.ref_start += self.FIXED_SHRINK_LEN
            alignment.query_end -= self.FIXED_SHRINK_LEN
        # end if
        return alignment
    # end def


    def _trim_start_primer(self, alignment, primer_num, left=True):
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
            alignment.query_start += primer_len_in_read
        # end if

        return alignment
    # end def

    def _trim_end_primer(self, alignment, primer_num, left=True):
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
            alignment.query_end -= primer_len_in_read
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
        self._REPORT_DELAY = max(
            1,
            round(self.NUM_READS_TOTAL * 0.01)
        )
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
        while num_done_reads.value >= next_report_num.value:
            next_report_num.value = next_report_num.value + self._REPORT_DELAY
        # end while
    # end def increment_next_report

    def print_status_bar(self):

        curr_num_done_reads = self.get_num_done_reads()

        bar_len = self._get_status_bar_len()
        ratio_done = curr_num_done_reads / self.NUM_READS_TOTAL
        progress_line_len = round(bar_len * ratio_done)

        print_arrow = progress_line_len != bar_len
        if print_arrow:
            arrow = '>'
            progress_line_len -= 1
        else:
            arrow = ''
        # end if

        sys.stdout.write(
            '\r{} - [{}{}{}] {}/{} ({}%)'.format(
                getwt(),
                '=' * progress_line_len,
                arrow,
                ' ' * (bar_len - progress_line_len),
                curr_num_done_reads,
                self.NUM_READS_TOTAL,
                round(ratio_done * 100)
            )
        )
        sys.stdout.flush()
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