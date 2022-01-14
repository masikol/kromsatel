
import gzip

from src.fastq import write_fastq_record


class UnpairedBinner:

    def __init__(self, unpaired_output):

        self.output = unpaired_output

        self.major_reads = list()
        self.minor_reads = list()
        self.uncertain_reads = list()
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

        outfpaths = (
            self.output.major_outfpath,
            self.output.minor_outfpath,
            self.output.uncertain_outfpath,
        )

        read_collections = (
            self.major_reads,
            self.minor_reads,
            self.uncertain_reads,
        )

        for outfpath, reads in zip(outfpaths, read_collections):
            if len(reads) != 0:
                self._append_to_outfile(reads, outfpath)
            # end if
        # end for
    # end def

    def _append_to_outfile(self, reads, outfpath):
        with gzip.open(outfpath, 'at') as outfile:
            for read in reads:
                write_fastq_record(read, outfile)
            # end for
        # end with
    # end def
# end class


class PairedBinner:

    def __init__(self, paired_output):

        self.output = paired_output

        self.major_forward_reads = list()
        self.major_reverse_reads = list()

        self.minor_forward_reads = list()
        self.minor_reverse_reads = list()

        self.uncertain_forward_reads = list()
        self.uncertain_reverse_reads = list()

        self.unpaired_forward_reads = list()
        self.unpaired_reverse_reads = list()
    # end def

    def add_major_pair(self, forward_read, reverse_read):
        self.major_forward_reads.append(forward_read)
        self.major_reverse_reads.append(reverse_read)
    # end def

    def add_minor_pair(self, forward_read, reverse_read):
        self.minor_forward_reads.append(forward_read)
        self.minor_reverse_reads.append(reverse_read)
    # end def

    def add_uncertain_pair(self, forward_read, reverse_read):
        self.uncertain_forward_reads.append(forward_read)
        self.uncertain_reverse_reads.append(reverse_read)
    # end def

    def add_forward_unpaired_read(self, read):
        self.unpaired_forward_reads.append(read)
    # end def

    def add_reverse_unpaired_read(self, read):
        self.unpaired_reverse_reads.append(read)
    # end def

    def write_binned_reads(self):

        outfpaths = (
            self.output.major_forward_outfpath,        self.output.major_reverse_outfpath,
            self.output.minor_forward_outfpath,        self.output.minor_reverse_outfpath,
            self.output.uncertain_forward_outfpath, self.output.uncertain_reverse_outfpath,
            self.output.unpaired_forward_outfpath,
            self.output.unpaired_reverse_outfpath,
        )

        read_collections = (
            self.major_forward_reads,        self.major_reverse_reads,
            self.minor_forward_reads,        self.minor_reverse_reads,
            self.uncertain_forward_reads, self.uncertain_reverse_reads,
            self.unpaired_forward_reads,
            self.unpaired_reverse_reads,
        )

        for outfpath, reads in zip(outfpaths, read_collections):
            if len(reads) != 0:
                self._append_to_outfile(reads, outfpath)
            # end if
        # end for
    # end def

    def _append_to_outfile(self, reads, outfpath):
        with gzip.open(outfpath, 'at') as outfile:
            for read in reads:
                write_fastq_record(read, outfile)
            # end for
        # end with
    # end def
# end class
