
import src.filesystem as fs
from src.fatal_errors import FatalError
from src.sequences import verify_sequence, get_non_iupac_chars


def read_fasta_sequence(file_path):

    with fs.open_file_may_by_gzipped(file_path, 'rt') as fasta_file:

        fasta_file.readline() # pass header
        sequence = ''
        line = fasta_file.readline().strip().upper()
        line_counter = 1

        while not (_is_header(line) or line == ''):
            line_counter += 1
            if not verify_sequence(line):
                non_iupac_chars = get_non_iupac_chars(line)
                error_msg = '\nError: a non-IUPAC character found' \
                            ' in line #{} of file `{}`.\n' \
                            'Bad characters are the following:\n  {}' \
                                .format(line_counter, file_path, ', '.join(non_iupac_chars))
                raise FatalError(error_msg)
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
