# -*- coding: utf-8 -*-

import re
import os
import sys
import getopt

from src.platform import platf_depend_exit


def handle_cl_args():
    # Function handles command line arguments.

    # Get arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvd:t:c:',
            ['help', 'version',
            'db=', 'threads=', 'chunk-size=',
            'am=', 'im=']
        )
    except getopt.GetoptError as opt_err:
        print( str(opt_err) )
        print('See help (`-h` option)')
        platf_depend_exit(2)
    # end try

    # Check and add fastq files:
    fq_fpaths = list()
    is_fastq = lambda f: not re.match(r'.*\.f(ast)?q(\.gz)?$', f) is None

    for arg in args:
        if not os.path.exists(arg):
            print('Error: file `{}` does not exist!'.format(arg))
            platf_depend_exit(1)
        elif not is_fastq(arg):
            print('Error: file `{}` does not look like a fastq file!'.format(arg))
            print('The script detects fastq files by extention and can process gzipped files.')
            print('Permitted extentions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`')
            platf_depend_exit(1)
        else:
            fq_fpaths.append(os.path.abspath(arg))
        # end if
    # end for

    # Check if there are any files to process
    if len(fq_fpaths) == 0:
        print('Please, specify input file(s).')
        print('Type `{} -h` to see help message.'.format(sys.argv[0]))
        platf_depend_exit(0)
    # end if

    args = {
        'fq_fpaths': fq_fpaths, # paths to input files
        'db_fpath': None,       # path to database
        'n_thr': 1,             # number of threads
        'chunk_size': 1000,     # fastq chunk
        'min_len_major': 100,   # minimum length of "major" read
        'min_len_minor': 25,    # minimum length of "minor" read
    }

    # Handle options
    for opt, arg in opts:
        # Get path to database
        if opt in ('-d', '--db'):
            if not os.path.exists( '{}.nhr'.format(arg) ):
                print('Error: database `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                args['db_fpath'] = os.path.abspath(arg)
            # end if

        # Handle chunk size
        elif opt in ('-c', '--chunk-size'):
            try:
                args['chunk_size'] = int(arg)
                if args['chunk_size'] <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: chunk_size (option `{}`) must be positive integer number.'\
                    .format(opt))
                platf_depend_exit(1)
            # end try

        # Handle number of threads
        elif opt in ('-t', '--threads'):
            try:
                args['n_thr'] = int(arg)
                if args['n_thr'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: number of threads must be integer number > 0!')
                print(' And here is your value: `{}`'.format(arg))
                platf_depend_exit(1)
            # end try
            if args['n_thr'] > len(os.sched_getaffinity(0)):
                print('''\nWarning! You have specified {} threads to use
      although {} are available.'''.format(args['n_thr'], len(os.sched_getaffinity(0))))
                args['n_thr'] = len(os.sched_getaffinity(0))
                print('Switched to {} threads.\n'.format(args['n_thr']))
            # end if

        # Handle min major read length
        elif opt == '--am':
            try:
                args['min_len_major'] = int(arg)
                if args['min_len_major'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid minimum major read length passed with option `{}`: {}'\
                    .format(opt, arg))
                print('It must be integer number > 0.')
                platf_depend_exit(1)
            # end try

        # Handle min minor read length
        elif opt == '--im':
            try:
                args['min_len_minor'] = int(arg)
                if args['min_len_minor'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid minimum minor read length passed with option `{}`: {}'\
                    .format(opt, arg))
                print('It must be integer number > 0.')
                platf_depend_exit(1)
            # end try
        # end if
    # end for

    if args['db_fpath'] is None:
        print('Error: option `-d` is mandatory.')
        print('Type {} -h to see help message.'.format(sys.argv[0]))
        platf_depend_exit(1)
    # end if

    return args
# end def handle_cl_args
