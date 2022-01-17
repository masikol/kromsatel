
import gzip

from src.fastq import write_fastq_record
from src.output import UnpairedOutput, PairedOutput


class UnpairedBinner:

    def __init__(self, outdir_path, output_prefix):

        self.output = UnpairedOutput(outdir_path, output_prefix)

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

    def add_major_read(self, read):
        self.major_reads.append(read)
    # end def

    def add_minor_read(self, read):
        self.minor_reads.append(read)
    # end def

    def add_uncertain_read(self, read):
        self.uncertain_reads.append(read)
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
# end class


class PairedBinner:

    def __init__(self, outdir_path, output_prefix):

        self.output = PairedOutput(outdir_path, output_prefix)

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

    def add_major_pair(self, frw_read, rvr_read):
        self.major_frw_reads.append(frw_read)
        self.major_rvr_reads.append(rvr_read)
    # end def

    def add_minor_pair(self, frw_read, rvr_read):
        self.minor_frw_reads.append(frw_read)
        self.minor_rvr_reads.append(rvr_read)
    # end def

    def add_uncertain_pair(self, frw_read, rvr_read):
        self.uncertain_frw_reads.append(frw_read)
        self.uncertain_rvr_reads.append(rvr_read)
    # end def

    def add_frw_unpaired_read(self, read):
        self.unpaired_frw_reads.append(read)
    # end def

    def add_rvr_unpaired_read(self, read):
        self.unpaired_rvr_reads.append(read)
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
# end class
