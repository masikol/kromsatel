
import src.fasta
import src.sequences
from src.printing import getwt
from src.alignment import Alignment
from src.fatal_errors import FatalError
from src.orientation import Orientation


class PrimerScheme:

    def __init__(self, kromsatel_args):
        self.primers_fpath = kromsatel_args.primers_fpath
        self.reference_fpath = kromsatel_args.reference_fpath
        self.primer_ext_len = kromsatel_args.primer_ext_len

        self.max_primer_len = 0
        self.primer_pairs = self._parse_primers()
    # end def

    def find_left_primer_by_coord(self, coord):
        for i, pair in enumerate(self.primer_pairs):
            if self.check_coord_within_primer(coord, i, Orientation.LEFT):
                return i
            # end if
        # end for
        return None
    # end def

    def find_right_primer_by_coord(self, coord):
        for i, pair in enumerate(self.primer_pairs):
            if self.check_coord_within_primer(coord, i, Orientation.RIGHT):
                return i
            # end if
        # end for
        return None
    # end def

    def check_coord_within_primer(self, coord, primer_pair_number, orientation):

        if primer_pair_number < 0 or primer_pair_number > len(self.primer_pairs)-1:
            return False
        # end if

        primer = self.get_primer(primer_pair_number, orientation)

        return primer.start <= coord <= primer.end
    # end def

    def get_primer(self, primer_pair_number, orientation):

        try:
            if orientation == Orientation.LEFT:
                primer = self.primer_pairs[primer_pair_number].left_primer
            else:
                primer = self.primer_pairs[primer_pair_number].right_primer
            # end if
        except TypeError:
            return None
        # end try

        return primer
    # end def

    def _parse_primers(self):

        sep = ','
        primer_pairs = list()

        n_lines = sum(1 for _ in open(self.primers_fpath, 'rt'))
        n_lines_is_even = n_lines % 2 == 0
        if not n_lines_is_even:
            error_msg = '\nError: Cannot parse primers from file `{}`.\n' \
                        'There are {} lines in this file.\n' \
                        'There must be even number of lines ' \
                        '(and therefore even number of primers), though.' \
                            .format(self.primers_fpath, self.primers_fpath)
            raise FatalError(error_msg)
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
                        Orientation.LEFT,
                        beg=find_start_pos
                    )
                    find_start_pos = left_start

                    right_start, right_end = self._find_primer_anneal_coords(
                        src.sequences.reverse_complement(right_primer_seq),
                        reference_seq,
                        Orientation.RIGHT,
                        beg=find_start_pos
                    )

                    primer_pairs.append(
                        PrimerPair(
                            Primer(left_start, left_end),
                            Primer(right_start, right_end),
                        )
                    )
                except ValueError as err:
                    error_msg = '\nError: cannot parse a line in file `{}`.\n{}' \
                        .format(self.primers_fpath, err)
                    raise FatalError(error_msg)
                # end try
            # end for
        # end with

        print('{} - Primers: found annealing coordinates'.format(getwt()))

        return primer_pairs
    # end def


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
    # end def

    def _parse_primer_from_csv_line(self, primer_line, sep=','):
        line_vals = primer_line.strip().split(sep)

        required_num_vals = 2
        if len(line_vals) < required_num_vals:
            error_msg = '\nError: not enough comma-separated columns.\n' \
                        '{} column(s) found, {} are required.\n' \
                        'The line: `{}`' \
                            .format(len(line_vals), required_num_vals, primer_line.strip())
            raise ValueError(error_msg)
        # end if

        primer_seq  = line_vals[1].upper()

        if not src.sequences.verify_sequence(primer_seq):
            error_msg = '\nError: a non-IUPAC character encountered' \
                        ' in the following primer sequence:\n' \
                        '  {}'.format(primer_seq)
            raise ValueError(error_msg)
        # end if

        return primer_seq
    # end def


    def _find_primer_anneal_coords(self, primer_seq, reference_seq, orientation, beg=0):

        start = reference_seq.find(primer_seq, beg)

        if start == -1:
            raise ValueError('Cannot find primer `{}` in the reference sequence'.format(primer_seq))
        # end if

        end = start + len(primer_seq) - 1

        if orientation == Orientation.LEFT:
            start = start - self.primer_ext_len
        else:
            end   =   end + self.primer_ext_len
        # end if

        return start, end
    # end def
# end class


class PrimerPair:
    def __init__(self, left_primer, right_primer):
        self.left_primer = left_primer
        self.right_primer = right_primer
    # end def

    def __repr__(self):
        return '[{}, {}]; [{}, {}];' \
            .format(
                self.left_primer.start,  self.left_primer.end,
                self.right_primer.start, self.right_primer.end
            )
    # end def
# end class


class Primer:
    def __init__(self, start, end):
        self.start = start # 0-based, left-closed
        self.end   = end   # 0-based, right-closed
    # end def
# end class
