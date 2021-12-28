
import os
import sys
import functools
import multiprocessing as mp

import src.fastq
import src.primers as prm
from src.progress import Progress
import src.synchronization as synchron
from src.printing import getwt, print_err
from src.alignment import parse_alignments, Alignment

from src.binning import PairedBinner
from src.binning import MAJOR, MINOR, NON_SPECIFIC



class ReadsCleaner:

    def __init__(self, args):
        self.args = args

        self.MIN_LEN = args['min_len']
        self.threads = args['n_thr']
        self.FIXED_CROP_LEN = args['fixed_crop_len']

        num_reads_total = self._count_reads()

        self.primer_scheme = prm.PrimerScheme(args)
        self.progress = Progress(num_reads_total)
    # end def __init__


    def clean_reads(self):

        # TODO
        # Unpaired mode
        if not self.args['paired_mode']:
            print_err('\nError: unpaired mode is not supported yet. :(')
            return
        # end if

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
    # end def clean_reads


    def _clean_reads_chunk_paired(self, reads_chunk):

        forward_chunk = reads_chunk[0]
        forward_alignments = parse_alignments(
            src.blast.blast_align(forward_chunk, self.args)
        )

        reverse_chunk = reads_chunk[1]
        reverse_alignments = parse_alignments(
            src.blast.blast_align(reverse_chunk, self.args)
        )

        binner = PairedBinner(self.args['output'])

        for forward_read, reverse_read in zip(*reads_chunk):

            forward_alignment = forward_alignments[forward_read['seq_id']]
            reverse_alignment = reverse_alignments[reverse_read['seq_id']]

            if forward_alignment is None or reverse_alignment is None:
                continue
            # end if

            forward_left = forward_alignment.align_strand_plus
            reverse_left = reverse_alignment.align_strand_plus

            improper_orientation = forward_left == reverse_left
            if improper_orientation:
                continue
            # end if

            classification = NON_SPECIFIC
            crop_forward_end, crop_reverse_end = True, True
            forward_survives, reverse_survives = False, False

            forward_start_primer_num, reverse_start_primer_num = None, None
            forward_end_primer_num,   reverse_end_primer_num   = None, None

            forward_start_coord = self._get_read_start_coord(forward_alignment)
            reverse_start_coord = self._get_read_start_coord(reverse_alignment)

            forward_start_primer_num = self._search_start_primer_bruteforce(
                forward_start_coord,
                forward_left
            )
            forward_start_primer_found = not (forward_start_primer_num is None)

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
                reverse_start_primer_num = self._search_start_primer_bruteforce(
                    reverse_start_coord,
                    forward_left
                )
            # end if

            if not forward_start_primer_num is None:
                reverse_end_coord = self._get_read_end_coord(reverse_alignment)
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
                forward_end_coord = self._get_read_end_coord(forward_alignment)
                crop_forward_end = self.primer_scheme.check_coord_within_primer(
                    forward_end_coord,
                    reverse_start_primer_num,
                    left=reverse_left
                )
                if crop_forward_end:
                    forward_end_primer_num = reverse_start_primer_num
                # end if
            # end if

            forward_alignment = self._trim_aligment(
                forward_alignment,
                forward_start_primer_num,
                forward_end_primer_num,
                crop_forward_end,
                forward_left
            )

            reverse_alignment = self._trim_aligment(
                reverse_alignment,
                reverse_start_primer_num,
                reverse_end_primer_num,
                crop_reverse_end,
                reverse_left
            )

            trimmed_forward_align_len = forward_alignment.get_align_len()
            trimmed_reverse_align_len = reverse_alignment.get_align_len()

            forward_survives = trimmed_forward_align_len >= self.MIN_LEN
            reverse_survives = trimmed_reverse_align_len >= self.MIN_LEN

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
                    binner.add_non_specific_pair(forward_read, reverse_read)
                # end if
            elif forward_survives and not reverse_survives:
                binner.add_forward_unpaired_read(forward_read)
            elif not forward_survives and reverse_survives:
                binner.add_reverse_unpaired_read(reverse_read)
            # end if
        # end for

        with synchron.output_lock:
            binner.write_binned_reads()
        # end with
        with synchron.status_update_lock:
            prev_next_value = self.progress.get_next_report_num()
            self.progress.increment_done(len(forward_chunk))
            self.progress.increment_next_report()
            if self.progress.get_next_report_num() != prev_next_value:
                with synchron.print_lock:
                    self.progress.print_status_bar()
                # end with
            # end if
        # end with
    # end def


    def _get_read_start_coord(self, alignment):
        if alignment.align_strand_plus:
            return alignment.ref_from
        else:
            return alignment.ref_to
        # end if
    # end def

    def _get_read_end_coord(self, alignment):
        if alignment.align_strand_plus:
            return alignment.ref_to
        else:
            return alignment.ref_from
        # end if
    # end def


    def _search_start_primer_bruteforce(self, start_coord, left=True):
        if left:
            start_primer_num = \
                self.primer_scheme.find_left_primer_by_coord(
                    start_coord
                )
        else:
            start_primer_num = \
                self.primer_scheme.find_right_primer_by_coord(
                    start_coord
                )
        # end if
        return start_primer_num
    # end def


    def _trim_aligment(self, alignment, start_primer_num, end_primer_num, crop_end=True, left=True):

        if not start_primer_num is None:
            alignment = self._trim_start_primer(alignment, start_primer_num, left)
        else:
            alignment = self._crop_start(alignment)
        # end if

        right = not left
        if not end_primer_num is None:
            alignment = self._trim_end_primer(alignment, end_primer_num, right)
        else:
            if crop_end:
                alignment = self._crop_end(alignment)
            # end if
        # end if

        return alignment
    # end def


    def _trim_read(self, read, alignment):
        new_start, new_end = alignment.query_from, alignment.query_to+1
        read['seq']  = read['seq'] [new_start : new_end]
        read['qual'] = read['qual'][new_start : new_end]
        return read
    # end def


    def _crop_start(self, alignment):
        if alignment.align_strand_plus:
            alignment.ref_from += self.FIXED_CROP_LEN
            alignment.query_from += self.FIXED_CROP_LEN
        else:
            alignment.ref_to -= self.FIXED_CROP_LEN
            alignment.query_from += self.FIXED_CROP_LEN
        # end if
        return alignment
    # end def

    def _crop_end(self, alignment):
        if alignment.align_strand_plus:
            alignment.ref_to -= self.FIXED_CROP_LEN
            alignment.query_to -= self.FIXED_CROP_LEN
        else:
            alignment.ref_from += self.FIXED_CROP_LEN
            alignment.query_to -= self.FIXED_CROP_LEN
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
            primer_len_in_read = primer.end - alignment.ref_from + 1
            alignment.ref_from += primer_len_in_read
            alignment.query_from += primer_len_in_read
        else:
            primer_len_in_read = alignment.ref_to - primer.start + 1
            alignment.ref_to -= primer_len_in_read
            alignment.query_from += primer_len_in_read
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
            primer_len_in_read = alignment.ref_to - primer.start + 1
            alignment.ref_to -= primer_len_in_read
            alignment.query_to -= primer_len_in_read
        else:
            primer_len_in_read = primer.end - alignment.ref_from + 1
            alignment.ref_from += primer_len_in_read
            alignment.query_to -= primer_len_in_read
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
