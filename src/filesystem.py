
import re
import os
import gzip
import shutil

from src.printing import print_err
from src.fatal_errors import FatalError


def open_file_may_by_gzipped(fpath, mode='rt'):
    file_is_gzipped = is_gzipped(fpath)

    if file_is_gzipped:
        open_func = gzip.open
    else:
        open_func = open
    # end if

    return open_func(fpath, mode)
# end def


def is_gzipped(fpath):
    if fpath.endswith('.gz'):
        try:
            with gzip.open(fpath) as _:
                pass
            # end with
        except gzip.BadGzipFile as err:
            error_msg = '\nError: bad gzip file: {}'.format(err)
            raise FatalError(error_msg)
        except OSError as err:
            error_msg = '\nError: cannot open file: {}'.format(err)
            raise FatalError(error_msg)
        # end try
    else:
        return False
    # end if

    return True
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
# end def


def init_file(fpath):
    try:
        with open(fpath, 'wt') as _:
            pass
        # end with
    except OSError as err:
        error_msg = '\nError: cannot initialize output file `{}`:\n {}' \
            .format(fpath, err)
        raise FatalError(error_msg)
    # end try
# end def


def create_dir(dirpath):
    if not os.path.exists(dirpath):
        try:
            os.makedirs(dirpath)
        except OSError as err:
            error_msg = '\nError: cannot create directory `{}`' \
                .format(dirpath)
            raise FatalError(error_msg)
        # end try
    # end if
# end def


def rm_file_warn_on_error(file_path):
    try:
        os.unlink(file_path)
    except OSError as oserr:
        print_err('\nWarning: Cannot remove file `{}`'.format(file_path))
        print_err( str(oserr) )
    # end try
# end def


def try_rm_directory(dirpath):
    try:
        shutil.rmtree(dirpath)
    except OSError as err:
        print_err('\nWarning: cannot remove temporary directory `{}`'.format(dirpath))
        print_err(str(err))
    # end try
# end def
