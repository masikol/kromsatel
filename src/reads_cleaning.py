
import sys
import multiprocessing as mp

import src.fastq
from src.printing import getwt
from src.progress import Progress
import src.synchronization as synchron
from src.binning import UnpairedBinner, PairedBinner
from src.classification import NanoporeReadsClassifier
from src.alignment import parse_alignments_illumina, parse_alignments_nanopore


class ReadsCleaner:

    def __init__(self, kromsatel_args):
        self.kromsatel_args = kromsatel_args
        self.threads_num = kromsatel_args.threads_num
    # end def

    def clean_reads(self):
        raise NotImplementedError
    # end def

    def _write_output(self, binner):
        with synchron.output_lock:
            binner.write_binned_reads()
        # end with
    # end def

    def _update_progress(self, increment):
        with synchron.status_update_lock:
            self.progress.increment_done(increment)
            self.progress.increment_next_report()
        # end with
    # end def

    def _print_progress(self):
        with synchron.print_lock:
            self.progress.print_status_bar()
        # end with
    # end def
# end class


class NanoporeReadsCleaner(ReadsCleaner):

    def __init__(self, kromsatel_args):
        super().__init__(kromsatel_args)
        self.classifier = NanoporeReadsClassifier(kromsatel_args)

        num_reads_total = \
            _count_reads_verbosely_unpaired(self.kromsatel_args.unpaired_read_fpath)
        self.progress = Progress(num_reads_total)
    # end def


    def clean_reads(self):

        reads_chunks = src.fastq.fastq_chunks_unpaired(
            fq_fpath=self.kromsatel_args.unpaired_read_fpath,
            chunk_size=self.kromsatel_args.chunk_size
        )

        self.progress.print_status_bar()

        self._clean_chunks(reads_chunks)

        self.progress.print_status_bar()
        print()
    # end def

    def _clean_chunks(self, reads_chunks):
        with mp.Pool(self.threads_num) as pool:
            task_iterator = pool.imap(
                self._clean_nanopore_chunk,
                reads_chunks,
                chunksize=1
            )
            for task in task_iterator:
                pass
            # end for
        # end with

        pool.close()
        pool.join()
    # end def

    def _clean_nanopore_chunk(self, reads_chunk):

        alignments = parse_alignments_nanopore(
            src.blast.blast_align(reads_chunk, self.kromsatel_args)
        )

        binner = UnpairedBinner(self.kromsatel_args.output)
        binner = self.classifier.fill_binner(reads_chunk, alignments, binner)

        self._write_output(binner)
        increment = len(reads_chunk)
        self._update_progress(increment)
        self._print_progress()
    # end def
# end class


def _count_reads_verbosely_unpaired(fastq_fpath):
    print('{} - Counting reads...'.format(getwt()))
    num_reads_total = src.fastq.count_reads(fastq_fpath)
    print('{} - {} reads.'.format(getwt(), num_reads_total))
    return num_reads_total
# end def


# class IlluminaPEReadsCleaner(ReadsCleaner):

#     def __init__(self, kromsatel_args):
#         super().__init__(kromsatel_args)
#         num_reads_total = self._count_reads()
#         self.progress = Progress(num_reads_total)
#     # end def


#     def clean_reads(self):

#         reads_chunks = src.fastq.fastq_chunks_paired(
#             forward_read_fpath=self.kromsatel_args.forward_read_fpath,
#             reverse_read_fpath=self.kromsatel_args.reverse_read_fpath,
#             chunk_size=self.kromsatel_args.chunk_size
#         )

#         self.progress.print_status_bar()

#         # Proceed
#         with mp.Pool(self.threads) as pool:
#             task_iterator = pool.imap(
#                 self._clean_paired_chunk,
#                 reads_chunks,
#                 chunksize=1
#             )
#             for task in task_iterator:
#                 pass
#             # end for
#         # end with

#         pool.close()
#         pool.join()

#         self.progress.print_status_bar()
#         print()
#     # end def


#     def _clean_paired_chunk(self, reads_chunk):

#         forward_chunk = reads_chunk[0]
#         forward_alignments = parse_alignments_illumina(
#             src.blast.blast_align(forward_chunk, self.kromsatel_args)
#         )

#         reverse_chunk = reads_chunk[1]
#         reverse_alignments = parse_alignments_illumina(
#             src.blast.blast_align(reverse_chunk, self.kromsatel_args)
#         )

#         binner = PairedBinner(self.kromsatel_args.output)

#         for forward_read, reverse_read in zip(*reads_chunk):

#             forward_alignment = forward_alignments[forward_read.header]
#             reverse_alignment = reverse_alignments[reverse_read.header]

#             if forward_alignment is None or reverse_alignment is None:
#                 continue
#             # end if

#             forward_orientation = self._get_orientation(forward_alignment)
#             reverse_orientation = self._get_orientation(reverse_alignment)

#             improper_orientation = forward_orientation == reverse_orientation
#             if improper_orientation:
#                 continue
#             # end if

#             classification = UNCERTAIN
#             crop_forward_end, crop_reverse_end = True, True
#             forward_survives, reverse_survives = False, False

#             forward_start_primer_num, reverse_start_primer_num = None, None
#             forward_end_primer_num,   reverse_end_primer_num   = None, None

#             forward_start_coord = self._get_read_start_coord(forward_alignment)
#             reverse_start_coord = self._get_read_start_coord(reverse_alignment)

#             forward_start_primer_num = self._search_primer_bruteforce(
#                 forward_start_coord,
#                 forward_orientation
#             )
#             forward_start_primer_found = not (forward_start_primer_num is None)

#             if forward_start_primer_found:

#                 reverse_start_major = \
#                     self.primer_scheme.check_coord_within_primer(
#                         reverse_start_coord,
#                         forward_start_primer_num,
#                         reverse_orientation
#                     )

#                 if reverse_start_major:
#                     classification = MAJOR
#                     reverse_start_primer_num = forward_start_primer_num
#                 else:
#                     minor_pair_primer_num = self._get_minor_primer_num(forward_start_primer_num, forward_orientation)
#                     reverse_start_minor = \
#                         self.primer_scheme.check_coord_within_primer(
#                             reverse_start_coord,
#                             minor_pair_primer_num,
#                             reverse_orientation
#                         )
#                     if reverse_start_minor:
#                         classification = MINOR
#                         reverse_start_primer_num = minor_pair_primer_num
#                     # end if
#                 # end if
#             else:
#                 reverse_start_primer_num = self._search_primer_bruteforce(
#                     reverse_start_coord,
#                     forward_orientation
#                 )
#             # end if

#             if not forward_start_primer_num is None:
#                 reverse_end_coord = self._get_read_end_coord(reverse_alignment)
#                 crop_reverse_end = self.primer_scheme.check_coord_within_primer(
#                     reverse_end_coord,
#                     forward_start_primer_num,
#                     forward_orientation
#                 )
#                 if crop_reverse_end:
#                     reverse_end_primer_num = forward_start_primer_num
#                 # end if
#             # end if
#             if not reverse_start_primer_num is None:
#                 forward_end_coord = self._get_read_end_coord(forward_alignment)
#                 crop_forward_end = self.primer_scheme.check_coord_within_primer(
#                     forward_end_coord,
#                     reverse_start_primer_num,
#                     reverse_orientation
#                 )
#                 if crop_forward_end:
#                     forward_end_primer_num = reverse_start_primer_num
#                 # end if
#             # end if

#             forward_alignment = self._trim_aligment(
#                 forward_alignment,
#                 forward_start_primer_num,
#                 forward_end_primer_num,
#                 forward_orientation,
#                 crop_forward_end
#             )

#             reverse_alignment = self._trim_aligment(
#                 reverse_alignment,
#                 reverse_start_primer_num,
#                 reverse_end_primer_num,
#                 reverse_orientation,
#                 crop_reverse_end
#             )

#             trimmed_forward_align_len = forward_alignment.get_align_len()
#             trimmed_reverse_align_len = reverse_alignment.get_align_len()

#             forward_survives = trimmed_forward_align_len >= self.MIN_LEN
#             reverse_survives = trimmed_reverse_align_len >= self.MIN_LEN

#             if forward_survives:
#                 forward_read = self._trim_read(forward_read, forward_alignment)
#             # end if
#             if reverse_survives:
#                 reverse_read = self._trim_read(reverse_read, reverse_alignment)
#             # end if

#             # Binning
#             if forward_survives and reverse_survives:
#                 if classification == MAJOR:
#                     binner.add_major_pair(forward_read, reverse_read)
#                 elif classification == MINOR:
#                     binner.add_minor_pair(forward_read, reverse_read)
#                 else:
#                     binner.add_uncertain_pair(forward_read, reverse_read)
#                 # end if
#             elif forward_survives and not reverse_survives:
#                 binner.add_forward_unpaired_read(forward_read)
#             elif not forward_survives and reverse_survives:
#                 binner.add_reverse_unpaired_read(reverse_read)
#             # end if
#         # end for

#         increment = len(forward_chunk)
#         self._write_output(binner)
#         self._update_progress(increment)
#         self._print_progress()
#     # end def

#     def _trim_aligment(self, alignment, start_primer_num, end_primer_num, orientation, crop_end=True):

#         if not start_primer_num is None:
#             alignment = self._trim_start_primer(alignment, start_primer_num, orientation)
#         else:
#             alignment = self._crop_start(alignment)
#         # end if

#         opposite_orientation = switch_orientation[orientation]
#         if not end_primer_num is None:
#             alignment = self._trim_end_primer(alignment, end_primer_num, opposite_orientation)
#         else:
#             if crop_end:
#                 alignment = self._crop_end(alignment)
#             # end if
#         # end if

#         return alignment
#     # end def

#     def _count_reads(self):
#         print('{} - Counting reads...'.format(getwt()))
#         num_reads_total = src.fastq.count_reads(self.kromsatel_args.forward_read_fpath)
#         print('{} - {} read pairs.'.format(getwt(), num_reads_total))
#         return num_reads_total
#     # end def
# # end class
