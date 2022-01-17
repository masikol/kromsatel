
import src.primers as prm
from src.orientation import Orientation
from src.trimming import NanoporeTrimmer, NanoporeTrimmingRule
from src.trimming import IlluminaPETrimmer, IlluminaPETrimmingRule
from src.trimming import ReadEndTrimmingRulePE
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

                alignment_overlaps = self._check_overlap(
                    alignment,
                    non_ovl_query_spans
                )

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

        read_start_coord = alignment.get_aln_start_coord()
        read_end_coord   = alignment.get_aln_end_coord()

        start_primer_num = self._search_primer_bruteforce(
            read_start_coord,
            read_orientation
        )

        start_primer_found = not (start_primer_num is None)

        classification_mark = self.UNCERTAIN

        if start_primer_found:
            if self._alignment_is_major(read_end_coord,
                                        start_primer_num,
                                        read_orientation):
                classification_mark = self.MAJOR
                end_primer_num = start_primer_num
            elif self._alignment_is_minor(read_end_coord,
                                          start_primer_num,
                                          read_orientation):
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

        trimming_rule = self._make_trimming_rule(
            start_primer_num,
            end_primer_num,
            read_orientation
        )

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


class IlluminaPEReadsClassifier(ReadsClassifier):

    def __init__(self, kromsatel_args):
        super().__init__(kromsatel_args)
        self.trimmer = IlluminaPETrimmer(kromsatel_args, self.primer_scheme)
    # end def

    def fill_binner(self, reads_chunk, alignments, binner):

        frw_alignments, rvr_alignments = alignments

        for frw_read, rvr_read in zip(*reads_chunk):

            frw_alignment = frw_alignments[frw_read.header]
            rvr_alignment = rvr_alignments[rvr_read.header]

            # TODO: process unpaired reads anyway
            if frw_alignment is None or rvr_alignment is None:
                continue
            # end if

            try:
                classification_mark, trimming_rules = \
                    self._classify_read_pair(frw_alignment, rvr_alignment)
            except ImproperOrientationError:
                continue
            # end try

            frw_alignment = \
                self.trimmer.trim_aligment(frw_alignment, trimming_rules[0])
            rvr_alignment = \
                self.trimmer.trim_aligment(rvr_alignment, trimming_rules[1])

            frw_trimmed_read = \
                self.trimmer.trim_read_to_fit_alignment(frw_read, frw_alignment)
            rvr_trimmed_read = \
                self.trimmer.trim_read_to_fit_alignment(rvr_read, rvr_alignment)

            self._add_read_pair_to_binner(
                frw_trimmed_read,
                rvr_trimmed_read,
                classification_mark,
                binner
            )
        # end for
    # end def

    def _classify_read_pair(self, frw_alignment, rvr_alignment):

        frw_orientation = get_read_orientation(frw_alignment)
        rvr_orientation = get_read_orientation(rvr_alignment)

        improper_orientation = (frw_orientation == rvr_orientation)
        if improper_orientation:
            raise ImproperOrientationError
        # end if

        frw_start_coord = frw_alignment.get_aln_start_coord()
        rvr_start_coord = rvr_alignment.get_aln_start_coord()

        frw_start_primer_num = self._search_primer_bruteforce(
            frw_start_coord,
            frw_orientation
        )
        frw_start_primer_found = not (frw_start_primer_num is None)

        classification_mark = self.UNCERTAIN

        if frw_start_primer_found:
            if self._alignments_are_major(rvr_start_coord,
                                          frw_start_primer_num,
                                          rvr_orientation):
                classification_mark = self.MAJOR
                rvr_start_primer_num = frw_start_primer_num
            elif self._alignments_are_minor(rvr_start_coord,
                                            frw_start_primer_num,
                                            frw_orientation,
                                            rvr_orientation):
                classification_mark = self.MINOR
                rvr_start_primer_num = _get_minor_primer_num(
                    frw_start_primer_num,
                    frw_orientation
                )
            # end if
        # end if

        if classification_mark == self.UNCERTAIN:
            rvr_start_primer_num = self._search_primer_bruteforce(
                rvr_start_coord,
                rvr_orientation
            )
        # end if

        frw_trimming_rule = self._make_trimming_rule(
            frw_alignment,
            frw_start_primer_num,
            frw_orientation,
            rvr_start_primer_num
        )
        rvr_trimming_rule = self._make_trimming_rule(
            rvr_alignment,
            rvr_start_primer_num,
            rvr_orientation,
            frw_start_primer_num
        )

        trimming_rules = (frw_trimming_rule, rvr_trimming_rule)

        return classification_mark, trimming_rules
    # end def

    def _alignments_are_major(self,
                              rvr_start_coord,
                              frw_start_primer_num,
                              rvr_orientation):
        alignments_are_major = \
            self.primer_scheme.check_coord_within_primer(
                rvr_start_coord,
                frw_start_primer_num,
                rvr_orientation
            )
        return alignments_are_major
    # end def

    def _alignments_are_minor(self,
                              rvr_start_coord,
                              frw_start_primer_num,
                              frw_orientation,
                              rvr_orientation):

        minor_pair_primer_num = _get_minor_primer_num(
            frw_start_primer_num,
            frw_orientation
        )

        alignments_are_minor = \
            self.primer_scheme.check_coord_within_primer(
                rvr_start_coord,
                minor_pair_primer_num,
                rvr_orientation
            )
        return alignments_are_minor
    # end def

    def _make_trimming_rule(self,
                            alignment,
                            start_primer_num,
                            read_orientation,
                            opposite_start_primer_num):

        start_primer = self.primer_scheme.get_primer(
            start_primer_num,
            read_orientation
        )

        end_trimming_rule = self._what_to_do_with_read_end(
            alignment,
            opposite_start_primer_num,
            read_orientation
        )

        return IlluminaPETrimmingRule(
            start_primer,
            read_orientation,
            end_trimming_rule
        )
    # end def

    def _what_to_do_with_read_end(self,
                                  alignment,
                                  opposite_start_primer_num,
                                  orientation):

        crop_end = True
        opposite_orientation = Orientation.invert(orientation)

        if not opposite_start_primer_num is None:
            end_coord = alignment.get_aln_end_coord()
            crop_end = not self.primer_scheme.check_coord_within_primer(
                end_coord,
                opposite_start_primer_num,
                opposite_orientation
            )
        # end if

        if not crop_end:
            end_primer_num = opposite_start_primer_num
            end_primer = self.primer_scheme.get_primer(
                end_primer_num,
                opposite_orientation
            )
        else:
            end_primer = None
        # end if

        return ReadEndTrimmingRulePE(end_primer, crop_end)
    # end def

    def _add_read_pair_to_binner(self,
                                 frw_read,
                                 rvr_read,
                                 classification_mark,
                                 binner):

        frw_long_enough = self._check_read_long_enough(frw_read)
        rvr_long_enough = self._check_read_long_enough(rvr_read)

        if frw_long_enough and rvr_long_enough:
            if classification_mark == self.MAJOR:
                binner.add_major_pair(frw_read, rvr_read)
            elif classification_mark == self.MINOR:
                binner.add_minor_pair(frw_read, rvr_read)
            else:
                binner.add_uncertain_pair(frw_read, rvr_read)
            # end if
        elif frw_long_enough and not rvr_long_enough:
            binner.add_frw_unpaired_read(frw_read)
        elif not frw_long_enough and rvr_long_enough:
            binner.add_rvr_unpaired_read(rvr_read)
        # end if
    # end def
# end class


class ImproperOrientationError(Exception):
    pass
# end class


def _get_minor_primer_num(primer_num, orientation):
    if orientation == Orientation.LEFT:
        minor_primer_num = primer_num - 1
    else:
        minor_primer_num = primer_num + 1
    # end if
    return minor_primer_num
# end def
