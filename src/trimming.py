
class TrimmingRule:

    def __init__(self, start_primer, orientation):
        self.start_primer = start_primer
        self.orientation = orientation
    # end def
# end class


class NanoporeTrimmingRule(TrimmingRule):

    def __init__(self, start_primer, end_primer, orientation):
        super().__init__(start_primer, orientation)
        self.end_primer = end_primer
    # end def
# end class


class IlluminaPETrimmingRule(TrimmingRule):

    def __init__(self, start_primer, orientation, read_end_trimming_rule):
        super().__init__(start_primer, orientation)
        self.read_end_trimming_rule = read_end_trimming_rule
    # end def
# end class


class ReadEndTrimmingRulePE:

    def __init__(self, end_primer, crop_end):
        self.end_primer = end_primer
        self.crop_end = crop_end
    # end def
# end class


class Trimmer:

    def __init__(self, kromsatel_args, primer_scheme):
        if kromsatel_args.fixed_crop_len == 'auto':
            self.FIXED_CROP_LEN = primer_scheme.max_primer_len
        else:
            self.FIXED_CROP_LEN = kromsatel_args.fixed_crop_len
        # end if
    # end def

    def trim_aligment(self, alignment, trimming_rule):
        raise NotImplementedError
    # end def

    def crop_start(self, alignment):
        if alignment.align_strand_plus:
            alignment.ref_from += self.FIXED_CROP_LEN
            alignment.query_from += self.FIXED_CROP_LEN
        else:
            alignment.ref_to -= self.FIXED_CROP_LEN
            alignment.query_from += self.FIXED_CROP_LEN
        # end if
        return alignment
    # end def

    def crop_end(self, alignment):
        if alignment.align_strand_plus:
            alignment.ref_to -= self.FIXED_CROP_LEN
            alignment.query_to -= self.FIXED_CROP_LEN
        else:
            alignment.ref_from += self.FIXED_CROP_LEN
            alignment.query_to -= self.FIXED_CROP_LEN
        # end if
        return alignment
    # end def

    def trim_start_primer(self, alignment, primer):

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

    def trim_end_primer(self, alignment, primer):

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

    def trim_read_to_fit_alignment(self, read, alignment):

        read_copy = read.get_copy()

        new_start, new_end = alignment.query_from, alignment.query_to+1

        read_copy.seq         = read_copy.seq[new_start : new_end]
        read_copy.quality_str = read_copy.quality_str[new_start : new_end]

        return read_copy
    # end def
# end class


class NanoporeTrimmer(Trimmer):

    def __init__(self, kromsatel_args, primer_scheme):
        super().__init__(kromsatel_args, primer_scheme)
    # end def

    def trim_aligment(self, alignment, trimming_rule):

        start_primer = trimming_rule.start_primer

        if not start_primer is None:
            alignment = self.trim_start_primer(alignment, start_primer)
        else:
            alignment = self.crop_start(alignment)
        # end if

        end_primer = trimming_rule.end_primer

        if not end_primer is None:
            alignment = self.trim_end_primer(alignment, end_primer)
        else:
            alignment = self.crop_end(alignment)
        # end if

        return alignment
    # end def
# end class


class IlluminaPETrimmer(Trimmer):

    def __init__(self, kromsatel_args, primer_scheme):
        super().__init__(kromsatel_args, primer_scheme)
    # end def

    def trim_aligment(self, alignment, trimming_rule):

        start_primer = trimming_rule.start_primer

        if not start_primer is None:
            alignment = self.trim_start_primer(alignment, start_primer)
        else:
            alignment = self.crop_start(alignment)
        # end if

        end_primer = trimming_rule.read_end_trimming_rule.end_primer

        # For Illumina, we do not necessarily need to crop end if end primer is not found --
        #   only if classification mark is UNCERTAIN
        if not end_primer is None:
            alignment = self.trim_end_primer(alignment, end_primer)
        else:
            if trimming_rule.read_end_trimming_rule.crop_end:
                alignment = self.crop_end(alignment)
            # end if
        # end if

        return alignment
    # end def
# end class
