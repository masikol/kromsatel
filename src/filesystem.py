
import re
import os
import gzip

from src.printing import print_err
from src.platform import platf_depend_exit


# For opening plain text and gzipped files
OPEN_FUNCS = (open, gzip.open)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode('utf-8').strip()  # format gzipped line
)


def rm_temp_file(file_path):
    try:
        os.unlink(file_path)
    except OSError as oserr:
        print_err('Warning: Cannot remove temporary file `{}`'.format(file_path))
        print_err( str(oserr) )
    # end try
# end def


def rm_fastq_extention(fpath):

    fastq_pattern = r'.+(\.f(ast)?q(\.gz)?)$'
    dirname = os.path.dirname(fpath)
    basename = os.path.basename(fpath)
    match_obj = re.match(fastq_pattern, basename)

    if not match_obj is None:
        extention = match_obj.group(1)
        new_basename = basename.replace(extention, '')
    else:
        new_basename = basename
    # end if

    return os.path.join(dirname, new_basename)
# end def rm_fastq_extention


def init_file(fpath):
    try:
        with open(fpath, 'wt') as _:
            pass
        # end with
    except OSError as err:
        print_err('\nError: cannot initialize output file `{}`'.format(fpath))
        print_err(str(err))
        platf_depend_exit(1)
    # end try
# end def init_file


def create_dir(dirpath):
    if not os.path.exists(dirpath):
        try:
            os.makedirs(dirpath)
        except OSError as err:
            print_err('\nError: cannot create directory `{}`'.format(dirpath))
            print_err(str(err))
            platf_depend_exit(1)
        # end try
    # end if
# end def create_dir
