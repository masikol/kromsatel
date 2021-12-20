


import src.blast
import src.output as out
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

    print(args)

    args['output'] = _configure_output(args)

    src.blast.create_reference_database(args)


    print('{} - Start.'.format(getwt()))

    reads_cleaner = rcl.ReadsCleaner(args)
    for x in reads_cleaner.primer_scheme.primer_pairs:
        print(x)

    return
    reads_cleaner.clean_reads()

    # print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
    # print('Output file: `{}`'.format(outfpath))
    # print('-------')

    print('{} - Completed.'.format(getwt()))
# end def main


def _configure_output(args):
    if args['paired_mode']:
        output = out.PairedOutput(args)
    else:
        output = out.UnpairedOutput(args)
    # end if
    return output
# end def _configure_output
