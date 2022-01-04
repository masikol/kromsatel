
import multiprocessing as mp

import src.fastq
import src.primers as prm
from src.printing import getwt
from src.progress import Progress
import src.synchronization as synchron
from src.alignment import parse_alignments_illumina, parse_alignments_nanopore, Alignment

from src.orientation import LEFT, RIGHT
from src.orientation import switch_orientation

from src.binning import UnpairedBinner, PairedBinner
from src.binning import MAJOR, MINOR, UNCERTAIN


class ReadsCleaner:

    def __init__(self, args):
        self.args = args

        self.MIN_LEN = args['min_len']
        self.threads = args['n_thr']

        self.primer_scheme = prm.PrimerScheme(args)

        if self.args['fixed_crop_len'] == 'auto':
            self.FIXED_CROP_LEN = self.primer_scheme.max_primer_len
        else:
            self.FIXED_CROP_LEN = args['fixed_crop_len']
        # end if
    # end def

    def clean_reads(self):
        raise NotImplementedError
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

    def _trim_start_primer(self, alignment, primer_num, orientation):
        if orientation == LEFT:
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

    def _trim_end_primer(self, alignment, primer_num, orientation):
        if orientation == LEFT:
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

    def _search_primer_bruteforce(self, coord, orientation):
        if orientation == LEFT:
            primer_num = \
                self.primer_scheme.find_left_primer_by_coord(
                    coord
                )
        else:
            primer_num = \
                self.primer_scheme.find_right_primer_by_coord(
                    coord
                )
        # end if
        return primer_num
    # end def

    def _trim_read(self, read, alignment):
        read_copy = read.copy()
        new_start, new_end = alignment.query_from, alignment.query_to+1
        read_copy['seq']  = read_copy['seq'] [new_start : new_end]
        read_copy['qual'] = read_copy['qual'][new_start : new_end]
        return read_copy
    # end def

    def _get_minor_primer_num(self, primer_num, orientation):
        if orientation == LEFT:
            minor_primer_num = primer_num - 1
        else:
            minor_primer_num = primer_num + 1
        # end if
        return minor_primer_num
    # end def

    def _get_orientation(self, alignment):
        if alignment.align_strand_plus:
            orientation = LEFT
        else:
            orientation = RIGHT
        # end if
        return orientation
    # end def

    def _write_output_and_print_progress(self, binner, increment):

        with synchron.output_lock:
            binner.write_binned_reads()
        # end with

        with synchron.status_update_lock:

            prev_next_value = self.progress.get_next_report_num()
            self.progress.increment_done(increment)
            self.progress.increment_next_report()

            if self.progress.get_next_report_num() != prev_next_value:
                with synchron.print_lock:
                    self.progress.print_status_bar()
                # end with
            # end if
        # end with
    # end def
# end class


class UnpairedReadsCleaner(ReadsCleaner):

    def __init__(self, args):
        super().__init__(args)
        num_reads_total = self._count_reads()
        self.progress = Progress(num_reads_total)
    # end def


    def clean_reads(self):

        reads_chunks = src.fastq.fastq_chunks_unpaired(
            fq_fpath=self.args['reads_unpaired'],
            chunk_size=self.args['chunk_size']
        )

        self.progress.print_status_bar()

        # Proceed
        with mp.Pool(self.threads) as pool:
            task_iterator = pool.imap(
                self._clean_unpaired_chunk,
                reads_chunks,
                chunksize=1
            )
            for task in task_iterator:
                pass
            # end for
        # end with

        pool.close()
        pool.join()

        self.progress.print_status_bar()
        print()
    # end def clean_reads


    def _clean_unpaired_chunk(self, reads_chunk):
        alignments = parse_alignments_nanopore(
            src.blast.blast_align(reads_chunk, self.args)
        )

        binner = UnpairedBinner(self.args['output'])

        for read in reads_chunk:

            read_alignments = alignments[read['seq_id']]

            if len(read_alignments) == 0:
                continue
            # end if

            non_ovl_query_spans = list()

            for alignment in read_alignments:

                alignment_overlaps = self._check_overlap(alignment, non_ovl_query_spans)

                if alignment_overlaps:
                    continue
                # end if
                non_ovl_query_spans.append(
                    (alignment.query_from, alignment.query_to,)
                )

                orientation = self._get_orientation(alignment)

                classification = UNCERTAIN

                read_survives = False

                start_primer_num = None
                end_primer_num = None

                read_start_coord = self._get_read_start_coord(alignment)
                read_end_coord   = self._get_read_end_coord(alignment)

                start_primer_num = self._search_primer_bruteforce(
                    read_start_coord,
                    orientation
                )
                start_primer_found = not (start_primer_num is None)

                if start_primer_found:

                    end_major = \
                        self.primer_scheme.check_coord_within_primer(
                            read_end_coord,
                            start_primer_num,
                            switch_orientation[orientation]
                        )

                    if end_major:
                        classification = MAJOR
                        end_primer_num = start_primer_num
                    else:
                        minor_pair_primer_num = self._get_minor_primer_num(start_primer_num, orientation)
                        end_minor = \
                            self.primer_scheme.check_coord_within_primer(
                                read_end_coord,
                                minor_pair_primer_num,
                                switch_orientation[orientation]
                            )
                        if end_minor:
                            classification = MINOR
                            end_primer_num = minor_pair_primer_num
                        # end if
                    # end if
                else:
                    end_primer_num = self._search_primer_bruteforce(
                        read_end_coord,
                        switch_orientation[orientation]
                    )
                # end if

                alignment = self._trim_aligment(
                    alignment,
                    start_primer_num,
                    end_primer_num,
                    orientation
                )

                trimmed_align_len = alignment.get_align_len()

                read_survives = trimmed_align_len >= self.MIN_LEN

                if read_survives:
                    read_to_write = self._trim_read(read, alignment)
                    read_to_write['seq_id'] = self._modify_read_name(read_to_write['seq_id'], alignment)

                    # Binning
                    if classification == MAJOR:
                        binner.add_major_read(read_to_write)
                    elif classification == MINOR:
                        binner.add_minor_read(read_to_write)
                    else:
                        binner.add_uncertain_read(read_to_write)
                    # end if
                # end if
            # end for

        # end for

        increment = len(reads_chunk)
        self._write_output_and_print_progress(binner, increment)
    # end def


    def _check_overlap(self, aligment, non_ovl_query_spans):
        for span in non_ovl_query_spans:
            span_from = span[0]
            span_to   = span[1]
            if span_from <= aligment.query_from <= span_to:
                return True
            # end if
            if span_from <= aligment.query_to   <= span_to:
                return True
            # end if
        # end for
        return False
    # end def


    def _count_reads(self):
        print('{} - Counting reads...'.format(getwt()))
        num_reads_total = src.fastq.count_reads(self.args['reads_unpaired'])
        print('{} - {} reads.'.format(getwt(), num_reads_total))
        return num_reads_total
    # end def


    def _trim_aligment(self, alignment, start_primer_num, end_primer_num, orientation):

        if not start_primer_num is None:
            alignment = self._trim_start_primer(alignment, start_primer_num, orientation)
        else:
            alignment = self._crop_start(alignment)
        # end if

        opposite_orientation = switch_orientation[orientation]
        if not end_primer_num is None:
            alignment = self._trim_end_primer(alignment, end_primer_num, opposite_orientation)
        else:
            alignment = self._crop_end(alignment)
        # end if

        return alignment
    # end def

    def _modify_read_name(self, read_name, alignment):
        identifier = read_name.partition(src.fastq.SPACE_HOLDER)[0]
        modified_identifier = '{}_{}-{}' \
            .format(identifier, alignment.query_from, alignment.query_to)
        return read_name.replace(identifier, modified_identifier)
    # end def
# end class


class PairedReadsCleaner(ReadsCleaner):

    def __init__(self, args):
        super().__init__(args)
        num_reads_total = self._count_reads()
        self.progress = Progress(num_reads_total)
    # end def


    def clean_reads(self):

        reads_chunks = src.fastq.fastq_chunks_paired(
            forward_read_fpath=self.args['reads_R1'],
            reverse_read_fpath=self.args['reads_R2'],
            chunk_size=self.args['chunk_size']
        )

        self.progress.print_status_bar()

        # Proceed
        with mp.Pool(self.threads) as pool:
            task_iterator = pool.imap(
                self._clean_paired_chunk,
                reads_chunks,
                chunksize=1
            )
            for task in task_iterator:
                pass
            # end for
        # end with

        pool.close()
        pool.join()

        self.progress.print_status_bar()
        print()
    # end def


    def _clean_paired_chunk(self, reads_chunk):

        forward_chunk = reads_chunk[0]
        forward_alignments = parse_alignments_illumina(
            src.blast.blast_align(forward_chunk, self.args)
        )

        reverse_chunk = reads_chunk[1]
        reverse_alignments = parse_alignments_illumina(
            src.blast.blast_align(reverse_chunk, self.args)
        )

        binner = PairedBinner(self.args['output'])

        for forward_read, reverse_read in zip(*reads_chunk):

            forward_alignment = forward_alignments[forward_read['seq_id']]
            reverse_alignment = reverse_alignments[reverse_read['seq_id']]

            if forward_alignment is None or reverse_alignment is None:
                continue
            # end if

            forward_orientation = self._get_orientation(forward_alignment)
            reverse_orientation = self._get_orientation(reverse_alignment)

            improper_orientation = forward_orientation == reverse_orientation
            if improper_orientation:
                continue
            # end if

            classification = UNCERTAIN
            crop_forward_end, crop_reverse_end = True, True
            forward_survives, reverse_survives = False, False

            forward_start_primer_num, reverse_start_primer_num = None, None
            forward_end_primer_num,   reverse_end_primer_num   = None, None

            forward_start_coord = self._get_read_start_coord(forward_alignment)
            reverse_start_coord = self._get_read_start_coord(reverse_alignment)

            forward_start_primer_num = self._search_primer_bruteforce(
                forward_start_coord,
                forward_orientation
            )
            forward_start_primer_found = not (forward_start_primer_num is None)

            if forward_start_primer_found:

                reverse_start_major = \
                    self.primer_scheme.check_coord_within_primer(
                        reverse_start_coord,
                        forward_start_primer_num,
                        reverse_orientation
                    )

                if reverse_start_major:
                    classification = MAJOR
                    reverse_start_primer_num = forward_start_primer_num
                else:
                    minor_pair_primer_num = self._get_minor_primer_num(forward_start_primer_num, forward_orientation)
                    reverse_start_minor = \
                        self.primer_scheme.check_coord_within_primer(
                            reverse_start_coord,
                            minor_pair_primer_num,
                            reverse_orientation
                        )
                    if reverse_start_minor:
                        classification = MINOR
                        reverse_start_primer_num = minor_pair_primer_num
                    # end if
                # end if
            else:
                reverse_start_primer_num = self._search_primer_bruteforce(
                    reverse_start_coord,
                    forward_orientation
                )
            # end if

            if not forward_start_primer_num is None:
                reverse_end_coord = self._get_read_end_coord(reverse_alignment)
                crop_reverse_end = self.primer_scheme.check_coord_within_primer(
                    reverse_end_coord,
                    forward_start_primer_num,
                    forward_orientation
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
                    reverse_orientation
                )
                if crop_forward_end:
                    forward_end_primer_num = reverse_start_primer_num
                # end if
            # end if

            forward_alignment = self._trim_aligment(
                forward_alignment,
                forward_start_primer_num,
                forward_end_primer_num,
                forward_orientation,
                crop_forward_end
            )

            reverse_alignment = self._trim_aligment(
                reverse_alignment,
                reverse_start_primer_num,
                reverse_end_primer_num,
                reverse_orientation,
                crop_reverse_end
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
                    binner.add_uncertain_pair(forward_read, reverse_read)
                # end if
            elif forward_survives and not reverse_survives:
                binner.add_forward_unpaired_read(forward_read)
            elif not forward_survives and reverse_survives:
                binner.add_reverse_unpaired_read(reverse_read)
            # end if
        # end for

        increment = len(forward_chunk)
        self._write_output_and_print_progress(binner, increment)
    # end def

    def _trim_aligment(self, alignment, start_primer_num, end_primer_num, orientation, crop_end=True):

        if not start_primer_num is None:
            alignment = self._trim_start_primer(alignment, start_primer_num, orientation)
        else:
            alignment = self._crop_start(alignment)
        # end if

        opposite_orientation = switch_orientation[orientation]
        if not end_primer_num is None:
            alignment = self._trim_end_primer(alignment, end_primer_num, opposite_orientation)
        else:
            if crop_end:
                alignment = self._crop_end(alignment)
            # end if
        # end if

        return alignment
    # end def

    def _count_reads(self):
        print('{} - Counting reads...'.format(getwt()))
        num_reads_total = src.fastq.count_reads(self.args['reads_R1'])
        print('{} - {} read pairs.'.format(getwt(), num_reads_total))
        return num_reads_total
    # end def
# end class
