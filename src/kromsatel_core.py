
import os
import multiprocessing as mp

import src.fastq
import src.filesystem as fs
from src.printing import getwt
from src.progress import Progress
import src.synchronization as synchron
from src.binning import SimpleUnpairedBinner, SimplePairedBinner, \
                        SplitUnpairedBinner, SplitPairedBinner
from src.alignment import parse_alignments_nanopore, parse_alignments_illumina
from src.reads_cleaning import NanoporeReadsCleaner, \
                               IlluminaPEReadsCleaner, \
                               IlluminaSEReadsCleaner


class KromsatelCore:

    def __init__(self, kromsatel_args):
        self.kromsatel_args = kromsatel_args
        self.threads_num = kromsatel_args.threads_num
    # end def

    def run(self):
        raise NotImplementedError
    # end def

    def _write_output(self):
        with synchron.output_lock:
            self.binner.write_binned_reads()
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


class LongReadKromsatelCore(KromsatelCore):

    # TODO: do not save reference to kromsatel args by creating some "BlastArguments" class
    def __init__(self, kromsatel_args):
        super().__init__(kromsatel_args)
        self.cleaner = NanoporeReadsCleaner(kromsatel_args)

        self.reads_fpath = self.kromsatel_args.long_read_fpath
        self.chunk_size = self.kromsatel_args.chunk_size

        num_reads_total = \
            _count_unpaired_reads_verbosely(self.reads_fpath)
        self.progress = Progress(num_reads_total)

        output_prefix = fs.rm_fastq_extention(
            os.path.basename(self.reads_fpath)
        )

        if kromsatel_args.split_output:
            self.binner = SplitUnpairedBinner(
                self.kromsatel_args.outdir_path,
                output_prefix,
                self.kromsatel_args.min_len
            )
        else:
            self.binner = SimpleUnpairedBinner(
                self.kromsatel_args.outdir_path,
                output_prefix,
                self.kromsatel_args.min_len
            )
        # end if
    # end def


    def run(self):

        reads_chunks = src.fastq.fastq_chunks_unpaired(
            fq_fpath=self.reads_fpath,
            chunk_size=self.chunk_size
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

        self.cleaner.fill_binner(reads_chunk, alignments, self.binner)

        self._write_output()
        increment = len(reads_chunk)
        self._update_progress(increment)
        self._print_progress()
    # end def
# end class


class IlluminaSEKromsatelCore(KromsatelCore):

    def __init__(self, kromsatel_args):
        super().__init__(kromsatel_args)
        self.cleaner = IlluminaSEReadsCleaner(kromsatel_args)

        self.reads_fpath = self.kromsatel_args.frw_read_fpath
        self.chunk_size = self.kromsatel_args.chunk_size

        num_reads_total = \
            _count_unpaired_reads_verbosely(self.reads_fpath)
        self.progress = Progress(num_reads_total)

        output_prefix = fs.rm_fastq_extention(
            os.path.basename(self.reads_fpath)
        )

        if kromsatel_args.split_output:
            self.binner = SplitUnpairedBinner(
                self.kromsatel_args.outdir_path,
                output_prefix,
                self.kromsatel_args.min_len
            )
        else:
            self.binner = SimpleUnpairedBinner(
                self.kromsatel_args.outdir_path,
                output_prefix,
                self.kromsatel_args.min_len
            )
        # end if
    # end def

    def run(self):

        reads_chunks = src.fastq.fastq_chunks_unpaired(
            fq_fpath=self.reads_fpath,
            chunk_size=self.chunk_size
        )

        self.progress.print_status_bar()

        self._clean_chunks(reads_chunks)

        self.progress.print_status_bar()
        print()
    # end def

    def _clean_chunks(self, reads_chunks):
        with mp.Pool(self.threads_num) as pool:
            task_iterator = pool.imap(
                self._clean_illumina_se_chunk,
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

    def _clean_illumina_se_chunk(self, reads_chunk):

        alignments = self._align_reads(reads_chunk)

        self.cleaner.fill_binner(reads_chunk, alignments, self.binner)

        self._write_output()
        increment = len(reads_chunk)
        self._update_progress(increment)
        self._print_progress()
    # end def

    def _align_reads(self, reads_chunk):
        alignments = parse_alignments_illumina(
            src.blast.blast_align(reads_chunk, self.kromsatel_args)
        )

        return alignments
    # end def
# end class


class IlluminaPEKromsatelCore(KromsatelCore):

    def __init__(self, kromsatel_args):
        super().__init__(kromsatel_args)
        self.cleaner = IlluminaPEReadsCleaner(kromsatel_args)

        self.frw_read_fpath = self.kromsatel_args.frw_read_fpath
        self.rvr_read_fpath = self.kromsatel_args.rvr_read_fpath
        self.chunk_size = self.kromsatel_args.chunk_size

        num_reads_total = \
            _count_paired_reads_verbosely(self.frw_read_fpath)
        self.progress = Progress(num_reads_total)

        output_prefix = fs.rm_fastq_extention(
            os.path.basename(self.frw_read_fpath)
        )

        if kromsatel_args.split_output:
            self.binner = SplitPairedBinner(
                self.kromsatel_args.outdir_path,
                output_prefix,
                self.kromsatel_args.min_len
            )
        else:
            self.binner = SimplePairedBinner(
                self.kromsatel_args.outdir_path,
                output_prefix,
                self.kromsatel_args.min_len
            )
        # end if
    # end def

    def run(self):

        reads_chunks = src.fastq.fastq_chunks_paired(
            frw_read_fpath=self.frw_read_fpath,
            rvr_read_fpath=self.rvr_read_fpath,
            chunk_size=self.chunk_size
        )

        self.progress.print_status_bar()

        self._clean_chunks(reads_chunks)

        self.progress.print_status_bar()
        print()
    # end def

    def _clean_chunks(self, reads_chunks):
        with mp.Pool(self.threads_num) as pool:
            task_iterator = pool.imap(
                self._clean_illumina_pe_chunk,
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

    def _clean_illumina_pe_chunk(self, reads_chunk):

        alignments = self._align_read_pairs(reads_chunk)

        self.cleaner.fill_binner(reads_chunk, alignments, self.binner)

        self._write_output()
        increment = len(reads_chunk[0])
        self._update_progress(increment)
        self._print_progress()
    # end def

    def _align_read_pairs(self, reads_chunk):
        frw_chunk = reads_chunk[0]
        frw_alignments = parse_alignments_illumina(
            src.blast.blast_align(frw_chunk, self.kromsatel_args)
        )

        rvr_chunk = reads_chunk[1]
        rvr_alignments = parse_alignments_illumina(
            src.blast.blast_align(rvr_chunk, self.kromsatel_args)
        )

        alignments = (frw_alignments, rvr_alignments)

        return alignments
    # end def
# end class


def _count_unpaired_reads_verbosely(fastq_fpath):
    print('{} - Counting reads...'.format(getwt()))
    num_reads_total = src.fastq.count_reads(fastq_fpath)
    print('{} - {} reads.'.format(getwt(), num_reads_total))
    return num_reads_total
# end def


def _count_paired_reads_verbosely(frw_fastq_fpath):
    print('{} - Counting reads...'.format(getwt()))
    num_reads_total = src.fastq.count_reads(frw_fastq_fpath)
    print('{} - {} read pairs.'.format(getwt(), num_reads_total))
    return num_reads_total
# end def
