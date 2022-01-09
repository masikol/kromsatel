
import re

import src.filesystem as fs
from src.printing import print_err
from src.platform import platf_depend_exit


def read_fasta_sequence(file_path):

    with fs.open_file_may_by_gzipped(file_path, 'rt') as fasta_file:

        fasta_file.readline() # pass header
        sequence = ''
        line = fasta_file.readline().strip().upper()
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
            line = fasta_file.readline().strip().upper()
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
