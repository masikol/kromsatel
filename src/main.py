# -*- encoding: utf-8 -*-

import src.blast
import src.shredding
import src.parse_args
import src.parse_primers
from src.printing import getwt


def main():
    # Handle command line arguments
    args = src.parse_args.handle_cl_args()

    # Check if blastn is installed
    src.blast.check_program('blastn')

    # Parse primers lengths if needed
    primers_lengths = None
    if not args['primers_fpath'] is None:
        primers_lengths = src.parse_primers.parse_primers_lengths(args['primers_fpath'])
    # end if

    print('{} - Start.'.format(getwt()))

    for fq_fpath in args['fq_fpaths']:
        print('{} - Processing file `{}`'.format(getwt(), fq_fpath))

        # Run cleaning
        outfpath = src.shredding.clean_and_shred(fq_fpath,
            args['db_fpath'],
            args['n_thr'],
            args['chunk_size'],
            args['min_len_major'],
            args['min_len_minor'],
            primers_lengths
        )

        print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
        print('Output file: `{}`'.format(outfpath))
        print('-------')
    # end for

    print('{} - Completed.'.format(getwt()))
# end def main
