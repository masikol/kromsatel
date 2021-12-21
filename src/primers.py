
import src.fasta
from src.platform import platf_depend_exit
from src.printing import print_err


_IUPAC_DICT = {
    'A': 'A',
    'T': 'T',
    'G': 'G',
    'C': 'C',
    'R': 'AG',
    'Y': 'CY',
    'W': 'AT',
    'S': 'GC',
    'K': 'GT',
    'M': 'AC',
    'H': 'ATC',
    'V': 'AGC',
    'B': 'TGC',
    'D': 'ATG',
    'N': 'ATGC',
}

_COMPLEMENT_DICT = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'R': 'Y',
    'Y': 'R',
    'W': 'W',
    'S': 'S',
    'K': 'M',
    'M': 'K',
    'H': 'D',
    'V': 'B',
    'B': 'V',
    'D': 'H',
    'N': 'N',
}


def _complement_base(base):
    return _COMPLEMENT_DICT[base]
# end def


def _reverse_complement(seq):

    return ''.join(map(_complement_base, seq))[::-1]
# end def


class PrimerScheme:

    def __init__(self, kromsatel_args):
        self.primers_fpath = kromsatel_args['primers_fpath']
        self.reference_fpath = kromsatel_args['reference_fpath']
        self.primer_pairs = self._parse_primers()
    # end def __init__

    # TODO
    # def find_left_primer(self, read):
    #     pass
    # # end def

    # TODO
    # def find_right_primer(self, read):
    #     pass
    # # end def

    def _parse_primers(self):

        sep = ','
        primer_pairs = list()

        n_lines = sum(1 for _ in open(self.primers_fpath, 'rt'))
        n_lines_is_even = n_lines % 2 == 0
        if not n_lines_is_even:
            print_err('Cannot parse primers: from file `{}`.'.format(self.primers_fpath))
            print_err('There are {} lines in this file.'.format(n_lines))
            print_err('There must be even number of lines \
    (and therefore even number of primers), though.')
            platf_depend_exit(1)
        # end if

        print('Parsing primers...')

        reference_seq = src.fasta.read_fasta_sequence(self.reference_fpath)
        find_start_pos = 0

        with open(self.primers_fpath, 'rt') as primers_file:
            for _ in range(n_lines // 2):
                try:

                    left_primer_seq, right_primer_seq = self._parse_primer_pair(primers_file, sep)

                    left_start, left_end = _find_primer_anneal_coords(
                        left_primer_seq,
                        reference_seq,
                        beg=find_start_pos
                    )
                    find_start_pos = left_start

                    right_start, right_end = _find_primer_anneal_coords(
                        _reverse_complement(right_primer_seq),
                        reference_seq,
                        beg=find_start_pos
                    )

                    primer_pairs.append(
                        PrimerPair(
                            Primer(left_start, left_end),
                            Primer(right_start, right_end),
                        )
                    )
                except ValueError as err:
                    print_err('Error: cannot parse a line in file `{}.`'.format(self.primers_fpath))
                    print_err(str(err))
                    platf_depend_exit(1)
                # end try
            # end for
        # end with

        print('Primers: annealing coordinates are found')

        return primer_pairs
    # end def parse_primers


    def _parse_primer_pair(self, primers_file, sep):


        left_primer_seq = self._parse_primer_from_csv_line(
            primers_file.readline(),
            sep
        )
        right_primer_seq = self._parse_primer_from_csv_line(
            primers_file.readline(),
            sep
        )

        return left_primer_seq, right_primer_seq
    # end def _parse_primer_pair

    def _parse_primer_from_csv_line(self, primer_line, sep=','):
        line_vals = primer_line.strip().split(sep)

        required_num_vals = 2
        if len(line_vals) != required_num_vals:
            raise ValueError(
                """Not enough comma-separated columns.
    {} column(s) found, {} are required'.
    The line: `{}`""".format(len(line_vals), required_num_vals, primer_line)
            )
        # end if

        primer_seq  = line_vals[1]

        return primer_seq
    # end def _parse_primer_from_csv_line

# end class PrimerScheme


def _find_primer_anneal_coords(primer_seq, reference_seq, beg=0):

    start = reference_seq.find(primer_seq, beg)

    if start == -1:
        raise ValueError('Cannot find primer `{}` in the reference sequence'.format(primer_seq))
    # end if

    end = start + len(primer_seq) - 1

    # Widen primer annealing interval by 1 bp in order not to
    #   misclassify minor alignments as non-specific ones due to single match
    #   occured by sheer chance.
    start = start - 1

    return start, end
# end def


class PrimerPair:
    def __init__(self, left_primer, right_primer):
        self.left_primer = left_primer
        self.right_primer = right_primer
    # end def __init__

    def __repr__(self):
        return '[{}, {}]; [{}, {}];' \
            .format(
                self.left_primer.start,  self.left_primer.end,
                self.right_primer.start, self.right_primer.end
            )
    # end def __repr__
# end class PrimerPair


class Primer:
    def __init__(self, start, end):
        self.start = start # 0-based, left-closed
        self.end   = end   # 0-based, right-closed
    # end def __init__
# end class Primer


def primer_match(read_seq, primer_seq, max_mismatch_num=2):

    num_mismatches = 0
    min_mismatch_to_discard = max_mismatch_num + 1
    local_iupac_dict = _IUPAC_DICT

    for pos, primer_base in enumerate(primer_seq):
        if not read_seq[pos] in local_iupac_dict[primer_base]:
            num_mismatches += 1
            if num_mismatches == min_mismatch_to_discard:
                return None
            # end if
        # end if
    # end for

    return pos
# end def primer_match
