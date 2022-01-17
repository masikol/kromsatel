
import os

import src.blast
import src.output as out
import src.filesystem as fs
import src.reads_cleaning as rcl
import src.parse_args
from src.printing import getwt, print_err
from src.platform import platformwise_exit

from src.fatal_errors import FatalError


def main():

    args = _parse_arguments()

    _check_blastplus_dependencies(args)

    # TODO: move output creation deeper
    output = _configure_output(args)
    args.set_output(output)

    db_fpath = src.blast.create_reference_database(args)
    args.set_database_path(db_fpath)

    print('{} - Start.'.format(getwt()))

    result_status = _clean_reads(args)

    _cleanup(args)

    if result_status == 0:
        print('\n{} - Completed.'.format(getwt()))
        print('  Output directory: `{}`'.format(args.outdir_path))
    else:
        print_err('\n{} - Completed with errors.'.format(getwt()))
    # end if
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


def _configure_output(kromsatel_args):
    if kromsatel_args.paired_mode:
        output = out.PairedOutput(kromsatel_args)
    else:
        output = out.UnpairedOutput(kromsatel_args)
    # end if
    return output
# end def


def _clean_reads(kromsatel_args):

    try:
        if kromsatel_args.paired_mode:
            reads_cleaner = rcl.IlluminaPEReadsCleaner(kromsatel_args)
        else:
            reads_cleaner = rcl.NanoporeReadsCleaner(kromsatel_args)
        # end if
        reads_cleaner.clean_reads()

    except FatalError as err:
        print_err(str(err))
        result_status = 1

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
