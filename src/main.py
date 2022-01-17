
import os

import src.blast
import src.parse_args
import src.filesystem as fs
import src.reads_cleaning as rcl
from src.printing import getwt, print_err
from src.platform import platformwise_exit
from src.fatal_errors import FatalError, InvalidFastqError


def main():

    args = _parse_arguments()

    _check_blastplus_dependencies(args)

    db_fpath = _create_database(args)
    args.set_database_path(db_fpath)

    print('{} - Start.'.format(getwt()))

    result_status = _clean_reads(args)

    _cleanup(args)

    if result_status == 0:
        print('\n{} - Completed.'.format(getwt()))
        print('  Output directory: `{}`'.format(args.outdir_path))
    else:
        print_err('\n\a{} - Completed with errors.'.format(getwt()))
    # end if

    return result_status
# end def


def _parse_arguments():
    try:
        args = src.parse_args.parse_args()
    except FatalError as err:
        print_err(str(err))
        platformwise_exit(1)
    # end try
    return args
# end def


def _check_blastplus_dependencies(kromsatel_args):
    blastplus_dependencies = src.blast.get_blastplus_dependencies(kromsatel_args)

    try:
        for ncbi_program in blastplus_dependencies:
            src.blast.check_program(ncbi_program)
        # end for
    except FatalError as err:
        print_err(str(err))
        platformwise_exit(1)
    # end try
# end def


def _create_database(kromsatel_args):
    try:
        db_fpath = src.blast.create_reference_database(kromsatel_args)
    except FatalError as err:
        print_err(str(err))
        platformwise_exit(1)
    # end try

    return db_fpath
# end def


def _clean_reads(kromsatel_args):

    result_status = 1

    try:
        if kromsatel_args.paired_mode:
            reads_cleaner = rcl.IlluminaPEReadsCleaner(kromsatel_args)
        else:
            reads_cleaner = rcl.NanoporeReadsCleaner(kromsatel_args)
        # end if
        reads_cleaner.clean_reads()

    except InvalidFastqError as err:
        print_err(err.msg_to_print)
        msg_to_log = '{}\n{}\n'.format(
            err.msg_to_print,
            err.msg_to_log_only
        )
        fs.log_to_file(msg_to_log, kromsatel_args.outdir_path)

    except FatalError as err:
        print_err(str(err))

    else:
        result_status = 0
    # end try

    return result_status
# end def


def _cleanup(kromsatel_args):
    fs.try_rm_directory(kromsatel_args.tmp_dir_path)

    fs.try_rm_directory(
        os.path.dirname(
            kromsatel_args.db_fpath
        )
    )
# end def
