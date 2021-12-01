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
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvd:t:c:p:o:',
            [
                'help', 'version',
                'db=', 'threads=', 'chunk-size=',
                'am=', 'im=', 'primers-to-rm=',
                'outdir='
            ]
        )
    except getopt.GetoptError as opt_err:
        print( str(opt_err) )
        print('See help (`-h` option)')
        platf_depend_exit(2)
    # end try

    # Check and add fastq files to list of files to process:
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

    kromsatel_args = {
        'fq_fpaths': fq_fpaths, # paths to input files
        'db_fpath': None,       # path to database
        'n_thr': 1,             # number of threads
        'chunk_size': 1000,     # fastq chunk
        'min_len_major': 100,   # minimum length of "major" read
        'min_len_minor': 25,    # minimum length of "minor" read
        'primers_fpath': None,    # path to file of primers to remove
        'outdir': os.path.join(os.getcwd(), 'kromsatel_output'), # output directory
    }

    # Handle options
    for opt, arg in opts:
        # Get path to database
        if opt in ('-d', '--db'):
            if not os.path.exists( '{}.nhr'.format(arg) ):
                print('Error: database `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                kromsatel_args['db_fpath'] = os.path.abspath(arg)
            # end if

        # Handle chunk size
        elif opt in ('-c', '--chunk-size'):
            try:
                kromsatel_args['chunk_size'] = int(arg)
                if kromsatel_args['chunk_size'] <= 0:
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
                kromsatel_args['n_thr'] = int(arg)
                if kromsatel_args['n_thr'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: number of threads must be integer number > 0!')
                print(' And here is your value: `{}`'.format(arg))
                platf_depend_exit(1)
            # end try
            if kromsatel_args['n_thr'] > len(os.sched_getaffinity(0)):
                print('''\nWarning! You have specified {} threads to use
      although {} are available.'''.format(kromsatel_args['n_thr'], len(os.sched_getaffinity(0))))
                kromsatel_args['n_thr'] = len(os.sched_getaffinity(0))
                print('Switched to {} threads.\n'.format(kromsatel_args['n_thr']))
            # end if

        # Handle min major read length
        elif opt == '--am':
            try:
                kromsatel_args['min_len_major'] = int(arg)
                if kromsatel_args['min_len_major'] < 1:
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
                kromsatel_args['min_len_minor'] = int(arg)
                if kromsatel_args['min_len_minor'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid minimum minor read length passed with option `{}`: {}'\
                    .format(opt, arg))
                print('It must be integer number > 0.')
                platf_depend_exit(1)
            # end try

        # Handle primers path
        elif opt in ('-p', '--primers-to-rm'):
            if not os.path.exists(arg):
                print('Error: file `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            # end if
            if not arg.endswith('.csv'):
                print('Error: primers file `{}` does not look like a csv file.'.format(arg))
                print("If the problem is just with file't extention, \
  please, change the extention to `.csv`.")
                platf_depend_exit(1)
            # end if
            kromsatel_args['primers_fpath'] = arg

        elif opt in ('-o', '--outdir'):
            kromsatel_args['outdir'] = os.path.abspath(arg)
            if not os.path.isdir(kromsatel_args['outdir']):
                try:
                    os.makedirs(kromsatel_args['outdir'])
                except OSError as err:
                    print('Error: cannot create output directory.')
                    print(str(err))
                    platf_depend_exit(1)
                # end try
            # end if
        # end if
    # end for

    if kromsatel_args['db_fpath'] is None:
        print('Error: option `-d` is mandatory.')
        print('Type {} -h to see help message.'.format(sys.argv[0]))
        platf_depend_exit(1)
    # end if

    primers_specified = not kromsatel_args['primers_fpath'] is None
    db_with_primers = 'with-primers' in kromsatel_args['db_fpath']

    if primers_specified and db_with_primers:
        print("Warning: you've specified option `-p`")
        print('  but your BLAST database `{}` looks like it contains primers.'\
            .format(kromsatel_args['db_fpath']))
        print("See https://github.com/masikol/kromsatel/blob/main/README.md#removing-primer-sequences for details")
        _ask_for_comfirmaton()
    elif not primers_specified and not db_with_primers:
        print("Warning: you haven't specified option `-p`")
        print("  but your BLAST database `{}` looks like it doesn't contain primers."\
            .format(kromsatel_args['db_fpath']))
        print("See https://github.com/masikol/kromsatel/blob/main/README.md#removing-primer-sequences for details")
        _ask_for_comfirmaton()
    # end if

    return kromsatel_args
# end def handle_cl_args


def _ask_for_comfirmaton():

    error = True
    while error:
        reply = input("""Press ENTER if you did it intentionally
  or enter `q` to quit: >> """)
        if reply == '':
            error = False
        elif reply.upper() == 'Q':
            print('Exitting...')
            sys.exit(0)
        else:
            print('Invalid reply: `{}`'.format(reply))
        # end if
    # end while
# end def _ask_for_comfirmaton
