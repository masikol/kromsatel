
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


class PrimerScheme:

    def __init__(self, primers_fpath):
        self.primers_fpath = primers_fpath
        self.primer_pairs = self._parse_primers()
    # end def __init__

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

        with open(self.primers_fpath, 'rt') as primers_file:
            for _ in range(n_lines // 2):
                try:
                    primer_pairs.append(
                        self._parse_primer_pair(primers_file, sep)
                    )
                except ValueError as err:
                    print_err('Error: cannot parse a line in file `{}.`'.format(self.primers_fpath))
                    print_err(str(err))
                    platf_depend_exit(1)
                # end try
            # end for
        # end with

        return primer_pairs
    # end def parse_primers


    def _parse_primer_pair(self, primers_file, sep):


        left_primer = self._parse_primer_from_csv_line(
            primers_file.readline(),
            sep
        )
        right_primer = self._parse_primer_from_csv_line(
            primers_file.readline(),
            sep
        )

        return PrimerPair(left_primer, right_primer)
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

        primer_name = line_vals[0]
        primer_seq  = line_vals[1]

        return Primer(primer_name, primer_seq)
    # end def _parse_primer_from_csv_line
# end class PrimerScheme


class PrimerPair:
    def __init__(self, left_primer, right_primer):
        self.left_primer = left_primer
        self.right_primer = right_primer
    # end def __init__

    def __repr__(self):
        return '{}:{}; {}:{};' \
            .format(
                self.left_primer.name, self.left_primer.seq,
                self.right_primer.name,
                self.right_primer.seq
            )
    # end def __repr__
# end class PrimerPair


class Primer:
    def __init__(self, name, seq):
        self.name = name
        self.seq  = seq
    # end def __init__
# end class Primer


def find_primer(read_seq, primer_seq, max_mismatch_num=2):

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

    return pos + 1
# end def find_primer
