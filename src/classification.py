
import src.primers as prm
from src.orientation import Orientation
from src.trimming import NanoporeTrimmer, NanoporeTrimmingRule
from src.orientation import get_read_orientation


class ReadsClassifier:

    MAJOR     = 0
    MINOR     = 1
    UNCERTAIN = 2

    def __init__(self, kromsatel_args):
        self.primer_scheme = prm.PrimerScheme(kromsatel_args)
        self.MIN_LEN = kromsatel_args.min_len
    # end def

    def fill_binner(self, reads_chunk, alignments, binner):
        raise NotImplementedError
    # end def

    def _search_primer_bruteforce(self, coord, orientation):
        if orientation == Orientation.LEFT:
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

    def _check_read_long_enough(self, read):
        return (len(read) > self.MIN_LEN)
    # end def
# enc class


class NanoporeReadsClassifier(ReadsClassifier):

    def __init__(self, kromsatel_args):
        super().__init__(kromsatel_args)
        self.trimmer = NanoporeTrimmer(kromsatel_args, self.primer_scheme)
    # end def

    def fill_binner(self, reads_chunk, alignments, binner):

        for read in reads_chunk:

            read_alignments = alignments[read.header]

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

                classification_mark, trimming_rule = \
                    self._classify_read(alignment)

                alignment = self.trimmer.trim_aligment(alignment, trimming_rule)

                trimmed_read = \
                    self.trimmer.trim_read_to_fit_alignment(read, alignment)

                if self._check_read_long_enough(trimmed_read):
                    trimmed_read = self._modify_read_name(trimmed_read, alignment)
                    self._add_read_to_binner(trimmed_read, classification_mark, binner)
                # end if
            # end for
        # end for

        return binner
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

    def _classify_read(self, alignment):

        read_orientation = get_read_orientation(alignment)

        read_start_coord = _get_read_start_coord(alignment)
        read_end_coord   = _get_read_end_coord(alignment)

        start_primer_num = self._search_primer_bruteforce(
            read_start_coord,
            read_orientation
        )

        start_primer_found = not (start_primer_num is None)

        classification_mark = self.UNCERTAIN

        if start_primer_found:
            if self._alignment_is_major(read_end_coord, start_primer_num, read_orientation):
                classification_mark = self.MAJOR
                end_primer_num = start_primer_num
            elif self._alignment_is_minor(read_end_coord, start_primer_num, read_orientation):
                classification_mark = self.MINOR
                end_primer_num = _get_minor_primer_num(start_primer_num, read_orientation)
            # end if
        # end if

        if classification_mark == self.UNCERTAIN:
            end_primer_num = self._search_primer_bruteforce(
                read_end_coord,
                Orientation.invert(read_orientation)
            )
        # end if

        trimming_rule = \
            self._make_trimming_rule(start_primer_num, end_primer_num, read_orientation)

        return classification_mark, trimming_rule
    # end def

    def _alignment_is_major(self, read_end_coord, start_primer_num, read_orientation):

        alignment_is_major = self.primer_scheme.check_coord_within_primer(
            read_end_coord,
            start_primer_num,
            Orientation.invert(read_orientation)
        )
        return alignment_is_major
    # end def

    def _alignment_is_minor(self, read_end_coord, start_primer_num, read_orientation):
        minor_pair_primer_num = _get_minor_primer_num(start_primer_num, read_orientation)
        alignment_is_minor = self.primer_scheme.check_coord_within_primer(
            read_end_coord,
            minor_pair_primer_num,
            Orientation.invert(read_orientation)
        )
        return alignment_is_minor
    # end def

    def _make_trimming_rule(self, start_primer_num, end_primer_num, read_orientation):
        start_primer = self.primer_scheme.get_primer(
            start_primer_num,
            read_orientation
        )
        end_primer = self.primer_scheme.get_primer(
            end_primer_num,
            Orientation.invert(read_orientation)
        )
        return NanoporeTrimmingRule(start_primer, end_primer, read_orientation)
    # end def

    def _modify_read_name(self, read, alignment):
        identifier = read.get_seqid()
        modified_identifier = \
            '{}_{}-{}'.format(identifier, alignment.query_from, alignment.query_to)

        read.header = read.header.replace(identifier, modified_identifier)
        return read
    # end def

    def _add_read_to_binner(self, read, classification_mark, binner):
        if classification_mark == self.MAJOR:
            binner.add_major_read(read)
        elif classification_mark == self.MINOR:
            binner.add_minor_read(read)
        else:
            binner.add_uncertain_read(read)
        # end if
    # end def
# end class


# TODO : move to Alignment class
def _get_read_start_coord(alignment):
    if alignment.align_strand_plus:
        return alignment.ref_from
    else:
        return alignment.ref_to
    # end if
# end def

# TODO : move to Alignment class
def _get_read_end_coord(alignment):
    if alignment.align_strand_plus:
        return alignment.ref_to
    else:
        return alignment.ref_from
    # end if
# end def


def _get_minor_primer_num(primer_num, orientation):
    if orientation == Orientation.LEFT:
        minor_primer_num = primer_num - 1
    else:
        minor_primer_num = primer_num + 1
    # end if
    return minor_primer_num
# end def
