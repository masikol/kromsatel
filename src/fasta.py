
import re

import src.filesystem as fs
from src.printing import print_err
from src.platform import platf_depend_exit


def read_fasta_sequence(file_path):

    is_gzipped = int(file_path.endswith('.gz'))
    open_func = fs.OPEN_FUNCS[is_gzipped]
    fmt_func = fs.FORMATTING_FUNCS[is_gzipped]

    with open_func(file_path) as fasta_file:

        fasta_file.readline() # pass header
        sequence = ''
        line = fmt_func(fasta_file.readline()).upper()
        line_counter = 1

        while not (_is_header(line) or line == ''):
            line_counter += 1
            if not _is_seq(line):
                print_err('A non-iupac character found in line #{} of file `{}`' \
                    .format(
                        line_counter, file_path
                    )
                )
                platf_depend_exit(1)
            # end if
            sequence += line
            line = fmt_func(fasta_file.readline()).upper()
        # end while
    # end with

    return sequence
# end def


def _is_header(string):
    return string.startswith('>')
# end def


def _is_seq(string):
    non_sequence_pattern = r'[^ATGCRYWSKMHVBDN]'
    is_valid_sequence = re.search(non_sequence_pattern, string) is None
    return is_valid_sequence
# end def
