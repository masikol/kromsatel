
import re


def parse_alignments(raw_alignments):

    alignments = dict()

    for raw_alignment_obj in raw_alignments:
        query_name = raw_alignment_obj['report']['results']['search']['query_title']
        alignments[query_name] = parse_single_alignment(raw_alignment_obj['report'])
    # end for

    return alignments
# end def


def parse_single_alignment(alignment_report):
    try:
        return Alignment(alignment_report)
    except NoAlignmentError:
        return None
    # end try
# end def


class Alignment:

    def __init__(self, blast_report):
        search_result = blast_report['results']['search']

        if len(search_result['hits']) == 0:
            raise NoAlignmentError()
        # end if

        hsp = search_result['hits'][0]['hsps'][0]

        self.query_name = search_result['query_title']
        self.align_len = hsp['align_len']

        self.query_start = hsp['query_from'] - 1 # 1-based to 0-based
        self.query_end   = hsp[ 'query_to' ] - 1 # 1-based to 0-based
        self.query_gap_locations = _find_gap_locations(hsp['qseq'])

        self.ref_start = hsp['hit_from'] - 1     # 1-based to 0-based
        self.ref_end   = hsp[ 'hit_to' ] - 1     # 1-based to 0-based
        self.ref_gap_locations = _find_gap_locations(hsp['hseq'])

        query_strand_plus = True if hsp['query_strand'].upper() == 'PLUS' else False
        ref_strand_plus   = True if hsp[ 'hit_strand' ].upper() == 'PLUS' else False
        self.align_strand_plus = query_strand_plus and ref_strand_plus

        # If strand is minus, swap reference alignment coordinates
        #   so that `ref_start` < `ref_end`
        if not self.align_strand_plus:
            self.ref_start, self.ref_end = self.ref_end, self.ref_start
        # end if
    # end def


    def __repr__(self):
        return '{}\nQ:[{}-{}];R:[{}-{}];({})\nQgaps:{}\nRgaps:{}' \
            .format(
                self.query_name,
                self.query_start,
                self.query_end,
                self.ref_start,
                self.ref_end,
                '+' if self.align_strand_plus else '-',
                self.query_gap_locations,
                self.ref_gap_locations
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
