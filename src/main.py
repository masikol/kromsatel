


# import src.blast
import src.output as out
# import src.shredding
import src.reads_cleaning as rcl
import src.parse_args
from src.printing import getwt
from src.platform import platf_depend_exit



def main():
    # Parse command line arguments
    args = src.parse_args.handle_cl_args()

    # Check if blastn is installed
    # src.blast.check_program('blastn')

    print(args)

    output = _configure_output(args)
    args['output'] = output

    print('{} - Start.'.format(getwt()))

    reads_cleaner = rcl.ReadsCleaner(args)
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

