
import gzip

from src.fastq import write_fastq_record


MAJOR        = 0
MINOR        = 1
NON_SPECIFIC = 2


class PairedBinner:

    def __init__(self, paired_output):

        self.output = paired_output

        self.major_forward_reads = list()
        self.major_reverse_reads = list()

        self.minor_forward_reads = list()
        self.minor_reverse_reads = list()

        self.non_specific_forward_reads = list()
        self.non_specific_reverse_reads = list()

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

    def add_non_specific_pair(self, forward_read, reverse_read):
        self.non_specific_forward_reads.append(forward_read)
        self.non_specific_reverse_reads.append(reverse_read)
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
            self.output.non_specific_forward_outfpath, self.output.non_specific_reverse_outfpath,
            self.output.unpaired_forward_outfpath,
            self.output.unpaired_reverse_outfpath,
        )

        read_collections = (
            self.major_forward_reads,        self.major_reverse_reads,
            self.minor_forward_reads,        self.minor_reverse_reads,
            self.non_specific_forward_reads, self.non_specific_reverse_reads,
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
