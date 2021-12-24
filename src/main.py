


import src.blast
import src.output as out
import src.filesystem as fs
# import src.shredding
import src.reads_cleaning as rcl
import src.parse_args
from src.printing import getwt
from src.platform import platf_depend_exit



def main():
    # Parse command line arguments
    args = src.parse_args.handle_cl_args()

    blastplus_dependencies = src.blast.get_blastplus_dependencies(args)

    for ncbi_program in blastplus_dependencies:
        src.blast.check_program(ncbi_program)
    # end for

    # print(args)

    args['output'] = _configure_output(args)

    args['db_fpath'] = src.blast.create_reference_database(args)

    print('{} - Start.'.format(getwt()))

    reads_cleaner = rcl.ReadsCleaner(args)
    # for i, x in enumerate(reads_cleaner.primer_scheme.primer_pairs):
    #     print(i, x)
    # return

    reads_cleaner.clean_reads()

    # print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
    # print('Output file: `{}`'.format(outfpath))
    # print('-------')

    _clean_tmp_files(args)

    print('\n{} - Completed.'.format(getwt()))
# end def


def _configure_output(kromsatel_args):
    if kromsatel_args['paired_mode']:
        output = out.PairedOutput(kromsatel_args)
    else:
        output = out.UnpairedOutput(kromsatel_args)
    # end if
    return output
# end def 


def _clean_tmp_files(kromsatel_args):
    fs.rm_tmp_dir(kromsatel_args['tmp_dir'])

    # TODO:
    # Uncomment. It is temporarily commented for debug purposes
    # fs.rm_tmp_dir(
    #     os.path.dirname(
    #         kromsatel_args['db_fpath']
    #     )
    # )
# end def
