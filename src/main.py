
import functools

import src.blast
import src.fastq
# import src.shredding
import src.parse_args
from src.printing import getwt
from src.primers import PrimerScheme
from src.platform import platf_depend_exit
from src.progress import Progress
from src.synchronizer import Synchronizer


def main():
    # Parse command line arguments
    args = src.parse_args.handle_cl_args()

    # Check if blastn is installed
    # src.blast.check_program('blastn')

    print(args)

    print('{} - Start.'.format(getwt()))

    outfpath = _clean_reads(args)

    print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
    print('Output file: `{}`'.format(outfpath))
    print('-------')

    print('{} - Completed.'.format(getwt()))
# end def main


def _clean_reads(args):

    # Parse primers
    args['primer_scheme'] = PrimerScheme(args['primers_fpath'])

    # Get path to output file.
    outfpath = src.filesystem.make_outfpath(fq_fpath, outdir)

    # Empty output file
    with open(outfpath, 'w') as _:
        pass
    # end with

    # Count reads and configure variables for printing status bar
    print('Counting reads...')
    num_reads_total = src.fastq.count_reads(fq_fpath)
    print('{} reads.'.format(num_reads_total))

    synchronizer = Synchronizer()

    progress = Progress(num_reads_total)
    progress.print_status_bar()

    fastq_chunks = _choose_fastq_chunks_func(args)

    # Proceed
    with mp.Pool(args['n_thr']) as pool:
        pool.starmap(_shredder, (
            (
                fq_chunk,
                args,
                progress,
                synchronizer
            )
            for fq_chunk in fastq_chunks())
        )
    # end with

    progress.print_status_bar()

    return outfpath
# end def _clean_reads


def _choose_fastq_chunks_func(args):
    if args['paired_mode']:
        fastq_chunks = functools.partial(
            src.fastq.fastq_chunks_paired,
            forward_read_fpath=args['reads_R1'],
            reverse_read_fpath=args['reads_R2'],
            chunk_size=args['chunk_size']
        )
    else:
        fastq_chunks = functools.partial(
            src.fastq.fastq_chunks_unpaired,
            fq_fpath=args['reads_unpaired'],
            chunk_size=args['chunk_size']
        )
    # end if
    return fastq_chunks
# end def _choose_fastq_chunks_func
