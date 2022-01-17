
import os

import src.blast
import src.output as out
import src.filesystem as fs
import src.reads_cleaning as rcl
import src.parse_args
from src.printing import getwt
from src.platform import platf_depend_exit



def main():
    args = src.parse_args.parse_args()

    blastplus_dependencies = src.blast.get_blastplus_dependencies(args)

    for ncbi_program in blastplus_dependencies:
        src.blast.check_program(ncbi_program)
    # end for

    # TODO: move output creation deeper
    output = _configure_output(args)
    args.set_output(output)

    db_fpath = src.blast.create_reference_database(args)
    args.set_database_path(db_fpath)

    print('{} - Start.'.format(getwt()))

    if args.paired_mode:
        reads_cleaner = rcl.IlluminaPEReadsCleaner(args)
    else:
        reads_cleaner = rcl.NanoporeReadsCleaner(args)
    # end if

    reads_cleaner.clean_reads()

    _clean_tmp_files(args)

    print('\n{} - Completed.'.format(getwt()))
    print('  Output directory: `{}`'.format(args.outdir_path))
# end def


def _configure_output(kromsatel_args):
    if kromsatel_args.paired_mode:
        output = out.PairedOutput(kromsatel_args)
    else:
        output = out.UnpairedOutput(kromsatel_args)
    # end if
    return output
# end def


def _clean_tmp_files(kromsatel_args):
    fs.try_rm_directory(kromsatel_args.tmp_dir_path)

    fs.try_rm_directory(
        os.path.dirname(
            kromsatel_args.db_fpath
        )
    )
# end def
