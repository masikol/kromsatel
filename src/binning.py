
import gzip

from src.fastq import write_fastq_record
from src.classification_marks import MAJOR, MINOR, UNCERTAIN
from src.output import SimpleUnpairedOutput, SimplePairedOutput, \
                       SplitUnpairedOutput, SplitPairedOutput


class Binner:

    def __init__(self, min_len):
        self.min_len = min_len
    # end def

    def write_binned_reads(self):
        raise NotImplementedError
    # end def

    def write_binned_reads(self):

        for outfpath, reads in zip(self.outfpaths, self.read_collections):
            if len(reads) != 0:
                self._append_to_outfile(reads, outfpath)
            # end if
        # end for

        self._clear()
    # end def

    def _append_to_outfile(self, reads, outfpath):
        with gzip.open(outfpath, 'at') as outfile:
            for read in reads:
                write_fastq_record(read, outfile)
            # end for
        # end with
    # end def

    def _clear(self):
        for collection in self.read_collections:
            collection.clear()
        # end for
    # end def

    def _check_read_long_enough(self, read):
        return (len(read) > self.min_len)
    # end def
# end class


class SimpleUnpairedBinner(Binner):

    def __init__(self, outdir_path, output_prefix, min_len):
        super().__init__(min_len)
        self.output = SimpleUnpairedOutput(outdir_path, output_prefix)
        self.output_reads = list()
    # end def

    def write_binned_reads(self):
        self._append_to_outfile(self.output_reads, self.output.outfpath)
        self._clear()
    # end def

    def _clear(self):
        self.output_reads.clear()
    # end def

    def add_read(self, read, classification_mark=None):
        if self._check_read_long_enough(read):
            self.output_reads.append(read)
        # end if
    # end def
# end class


class SimplePairedBinner(Binner):

    def __init__(self, outdir_path, output_prefix, min_len):

        super().__init__(min_len)
        self.output = SimplePairedOutput(outdir_path, output_prefix)

        self.frw_reads = list()
        self.rvr_reads = list()

        self.unpaired_frw_reads = list()
        self.unpaired_rvr_reads = list()

        self.outfpaths = (
            self.output.frw_outfpath,
            self.output.rvr_outfpath,
            self.output.frw_unpaired_outfpath,
            self.output.rvr_unpaired_outfpath,
        )

        self.read_collections = (
            self.frw_reads,
            self.rvr_reads,
            self.unpaired_frw_reads,
            self.unpaired_rvr_reads,
        )
    # end def

    def add_read_pair(self, frw_read, rvr_read, classification_mark=None):

        frw_long_enough = self._check_read_long_enough(frw_read)
        rvr_long_enough = self._check_read_long_enough(rvr_read)

        if frw_long_enough and rvr_long_enough:
            self._add_normal_read_pair(frw_read, rvr_read)
        elif frw_long_enough:
            self.unpaired_frw_reads.append(frw_read)
        elif rvr_long_enough:
            self.unpaired_rvr_reads.append(rvr_read)
        # end if
    # end def

    def _add_normal_read_pair(self, frw_read, rvr_read):
        self.frw_reads.append(frw_read)
        self.rvr_reads.append(rvr_read)
    # end def
# end class


class SplitUnpairedBinner(Binner):

    def __init__(self, outdir_path, output_prefix, min_len):

        super().__init__(min_len)
        self.output = SplitUnpairedOutput(outdir_path, output_prefix)

        self.major_reads = list()
        self.minor_reads = list()
        self.uncertain_reads = list()

        self.outfpaths = (
            self.output.major_outfpath,
            self.output.minor_outfpath,
            self.output.uncertain_outfpath,
        )

        self.read_collections = (
            self.major_reads,
            self.minor_reads,
            self.uncertain_reads,
        )
    # end def

    def add_read(self, read, classification_mark):

        if not self._check_read_long_enough(read):
            return
        # end if

        if classification_mark == MAJOR:
            self._add_major_read(read)
        elif classification_mark == MINOR:
            self._add_minor_read(read)
        else:
            self._add_uncertain_read(read)
        # end if
    # end def

    def _add_major_read(self, read):
        self.major_reads.append(read)
    # end def

    def _add_minor_read(self, read):
        self.minor_reads.append(read)
    # end def

    def _add_uncertain_read(self, read):
        self.uncertain_reads.append(read)
    # end def
# end class


class SplitPairedBinner(Binner):

    def __init__(self, outdir_path, output_prefix, min_len):

        super().__init__(min_len)
        self.output = SplitPairedOutput(outdir_path, output_prefix)

        self.major_frw_reads = list()
        self.major_rvr_reads = list()

        self.minor_frw_reads = list()
        self.minor_rvr_reads = list()

        self.uncertain_frw_reads = list()
        self.uncertain_rvr_reads = list()

        self.unpaired_frw_reads = list()
        self.unpaired_rvr_reads = list()

        self.outfpaths = (
            self.output.major_frw_outfpath,     self.output.major_rvr_outfpath,
            self.output.minor_frw_outfpath,     self.output.minor_rvr_outfpath,
            self.output.uncertain_frw_outfpath, self.output.uncertain_rvr_outfpath,
            self.output.unpaired_frw_outfpath,
            self.output.unpaired_rvr_outfpath,
        )

        self.read_collections = (
            self.major_frw_reads,     self.major_rvr_reads,
            self.minor_frw_reads,     self.minor_rvr_reads,
            self.uncertain_frw_reads, self.uncertain_rvr_reads,
            self.unpaired_frw_reads,
            self.unpaired_rvr_reads,
        )
    # end def

    def add_read_pair(self, frw_read, rvr_read, classification_mark):

        frw_long_enough = self._check_read_long_enough(frw_read)
        rvr_long_enough = self._check_read_long_enough(rvr_read)

        if frw_long_enough and rvr_long_enough:
            if classification_mark == MAJOR:
                self._add_major_pair(frw_read, rvr_read)
            elif classification_mark == MINOR:
                self._add_minor_pair(frw_read, rvr_read)
            else:
                self._add_uncertain_pair(frw_read, rvr_read)
            # end if
        elif frw_long_enough:
            self.unpaired_frw_reads.append(frw_read)
        elif rvr_long_enough:
            self.unpaired_rvr_reads.append(rvr_read)
        # end if
    # end def

    def _add_major_pair(self, frw_read, rvr_read):
        self.major_frw_reads.append(frw_read)
        self.major_rvr_reads.append(rvr_read)
    # end def

    def _add_minor_pair(self, frw_read, rvr_read):
        self.minor_frw_reads.append(frw_read)
        self.minor_rvr_reads.append(rvr_read)
    # end def

    def _add_uncertain_pair(self, frw_read, rvr_read):
        self.uncertain_frw_reads.append(frw_read)
        self.uncertain_rvr_reads.append(rvr_read)
    # end def
# end class
