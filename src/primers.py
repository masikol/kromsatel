
import src.fasta
from src.printing import print_err, getwt
from src.alignment import Alignment
from src.platform import platf_depend_exit


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


def _reverse_complement(seq):
    return ''.join(map(_complement_base, seq))[::-1]
# end def


def _complement_base(base):
    return _COMPLEMENT_DICT[base]
# end def


class PrimerScheme:

    def __init__(self, kromsatel_args):
        self.primers_fpath = kromsatel_args['primers_fpath']
        self.reference_fpath = kromsatel_args['reference_fpath']
        self.primer_ext_len = kromsatel_args['primer_ext_len']

        self.max_primer_len = 0
        self.primer_pairs = self._parse_primers()
    # end def __init__

    def find_left_primer_by_coord(self, coord):
        for i, pair in enumerate(self.primer_pairs):
            if self.check_coord_within_primer(coord, i, left=True):
                return i
            # end if
        # end for
        return None
    # end def

    def find_right_primer_by_coord(self, coord):
        for i, pair in enumerate(self.primer_pairs):
            if self.check_coord_within_primer(coord, i, left=False):
                return i
            # end if
        # end for
        return None
    # end def

    def check_coord_within_primer(self, coord, primer_pair_number, left=True):

        if primer_pair_number < 0 or primer_pair_number > len(self.primer_pairs)-1:
            return False
        # end if

        if left:
            primer = self.primer_pairs[primer_pair_number].left_primer
        else:
            primer = self.primer_pairs[primer_pair_number].right_primer
        # end if
        return primer.start <= coord <= primer.end
    # end def

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

        print('{} - Parsing primers...'.format(getwt()))

        reference_seq = src.fasta.read_fasta_sequence(self.reference_fpath)
        find_start_pos = 0

        with open(self.primers_fpath, 'rt') as primers_file:
            for _ in range(n_lines // 2):
                try:

                    left_primer_seq, right_primer_seq = self._parse_primer_pair(primers_file, sep)

                    self.max_primer_len = max(
                        self.max_primer_len,
                        len(left_primer_seq),
                        len(right_primer_seq)
                    )

                    left_start, left_end = self._find_primer_anneal_coords(
                        left_primer_seq,
                        reference_seq,
                        left=True,
                        beg=find_start_pos
                    )
                    find_start_pos = left_start

                    right_start, right_end = self._find_primer_anneal_coords(
                        _reverse_complement(right_primer_seq),
                        reference_seq,
                        left=False,
                        beg=find_start_pos
                    )

                    primer_pairs.append(
                        PrimerPair(
                            Primer(left_start, left_end),
                            Primer(right_start, right_end),
                        )
                    )
                except ValueError as err:
                    print_err('Error: cannot parse a line in file `{}`.'.format(self.primers_fpath))
                    print_err(str(err))
                    platf_depend_exit(1)
                # end try
            # end for
        # end with

        print('{} - Primers: found annealing coordinates'.format(getwt()))

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


    def _find_primer_anneal_coords(self, primer_seq, reference_seq, left=True, beg=0):

        start = reference_seq.find(primer_seq, beg)

        if start == -1:
            raise ValueError('Cannot find primer `{}` in the reference sequence'.format(primer_seq))
        # end if

        end = start + len(primer_seq) - 1

        if left:
            start = start - self.primer_ext_len
        else:
            end   =   end + self.primer_ext_len
        # end if

        return start, end
    # end def
# end class PrimerScheme




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
