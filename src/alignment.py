
import re


def parse_alignments_illumina(raw_alignments):

    alignments = dict()

    for raw_alignment_obj in raw_alignments:
        query_name = raw_alignment_obj['report']['results']['search']['query_title']
        alignments[query_name] = parse_single_hsp(raw_alignment_obj['report'])
    # end for

    return alignments
# end def


def parse_alignments_nanopore(raw_alignments):

    alignments = dict()

    for raw_alignment_obj in raw_alignments:
        query_name = raw_alignment_obj['report']['results']['search']['query_title']
        alignments[query_name] = parse_hsps(raw_alignment_obj['report'])
    # end for

    return alignments
# end def


def parse_single_hsp(alignment_report):

    search_result = alignment_report['results']['search']

    if len(search_result['hits']) == 0:
        return None
    # end if

    hsps = search_result['hits'][0]['hsps']

    first_hsp = hsps[0]
    return Alignment(first_hsp)
# end def


def parse_hsps(alignment_report):

    alignments = list()

    search_result = alignment_report['results']['search']

    if len(search_result['hits']) == 0:
        return alignments
    # end if

    hsps = search_result['hits'][0]['hsps']

    for hsp in hsps:
        alignments.append(Alignment(hsp))
    # end for

    return alignments
# end def


class Alignment:

    def __init__(self, hsp):

        self.query_from = hsp['query_from'] - 1 # 1-based to 0-based
        self.query_to   = hsp[ 'query_to' ] - 1 # 1-based to 0-based

        self.ref_from   = hsp[ 'hit_from' ] - 1 # 1-based to 0-based
        self.ref_to     = hsp[  'hit_to'  ] - 1 # 1-based to 0-based

        # TODO: Temporarily disabled
        # self.query_gap_locations = _find_gap_locations(hsp['qseq'])
        # self.ref_gap_locations = _find_gap_locations(hsp['hseq'])

        query_strand_plus = True if hsp['query_strand'].upper() == 'PLUS' else False
        ref_strand_plus   = True if hsp[ 'hit_strand' ].upper() == 'PLUS' else False
        self.align_strand_plus = query_strand_plus and ref_strand_plus

        # If strand is minus, swap reference alignment coordinates
        #   so that `ref_from` < `ref_to`
        if not self.align_strand_plus:
            self.ref_from, self.ref_to = self.ref_to, self.ref_from
        # end if
    # end def

    def get_align_len(self):
        return max(0, self.query_to - self.query_from)
    # end def

    def __repr__(self):
        # return '{}\nQ:[{}-{}];R:[{}-{}];({})\nQgaps:{}\nRgaps:{}' \
        return 'Q:[{}-{}];R:[{}-{}];({})' \
            .format(
                self.query_from,
                self.query_to,
                self.ref_from,
                self.ref_to,
                '+' if self.align_strand_plus else '-'#,
                # self.query_gap_locations,
                # self.ref_gap_locations
            )
    # end def
# end class


class NoAlignmentError(Exception):
    pass
# end class


def _find_gap_locations(sequence, gap='-'):
    gap_locations = tuple(
        (match.start() for match in re.finditer(gap, sequence))
    )
    return gap_locations
# end def
