# -*- coding: utf-8 -*-

import re
import os
import sys
import getopt

from src.printing import print_err
from src.platform import platf_depend_exit


def handle_cl_args():
    # Function handles command line arguments.

    # Get arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hv1:2:u:d:t:c:p:o:',
            [
                'help', 'version',
                'reads-R1=', 'reads-R2=', 'reads-unpaired=',
                'db=', 'threads=', 'chunk-size=',
                'am=', 'im=', 'primers=',
                'outdir='
            ]
        )
    except getopt.GetoptError as opt_err:
        print_err( str(opt_err) )
        print_err('See help (`-h` option)')
        platf_depend_exit(2)
    # end try

    # # Check and add fastq files to list of files to process:
    # fq_fpaths = list()
    # is_fastq = lambda f: not re.match(r'.*\.f(ast)?q(\.gz)?$', f) is None

    # for arg in args:
    #     if not os.path.exists(arg):
    #         print_err('Error: file `{}` does not exist!'.format(arg))
    #         platf_depend_exit(1)
    #     elif not is_fastq(arg):
    #         print_err('Error: file `{}` does not look like a fastq file!'.format(arg))
    #         print_err('The script detects fastq files by extention and can process gzipped files.')
    #         print_err('Permitted extentions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`')
    #         platf_depend_exit(1)
    #     else:
    #         fq_fpaths.append(os.path.abspath(arg))
    #     # end if
    # # end for

    # # Check if there are any files to process
    # if len(fq_fpaths) == 0:
    #     print_err('Please, specify input file(s).')
    #     print_err('Type `{} -h` to see help message.'.format(sys.argv[0]))
    #     platf_depend_exit(0)
    # # end if

    kromsatel_args = {
        'reads_R1': list(),
        'reads_R2': list(),
        'reads_unpaired': list(),
        'db_fpath': None,         # path to database
        'primers_fpath': None,    # path to file of primers to remove
        'n_thr': 1,               # number of threads
        'chunk_size': 1000,       # fastq chunk
        'min_len_major': 100,     # minimum length of "major" read
        'min_len_minor': 25,      # minimum length of "minor" read
        'outdir': os.path.join(os.getcwd(), 'kromsatel_output'), # output directory
        'paired_mode': False
    }

    # Handle options
    for opt, arg in opts:
        # Get path to database
        if opt in ('-1', '--reads-R1'):
            fpath = os.path.abspath(arg)
            try:
                _check_fastq_file(fpath)
            except ValueError as err:
                print_err(str(err))
                platf_depend_exit(1)
            else:
                kromsatel_args['reads_R1'].append(fpath)
            # end try

        elif opt in ('-2', '--reads-R2'):
            fpath = os.path.abspath(arg)
            try:
                _check_fastq_file(fpath)
            except ValueError as err:
                print_err(str(err))
                platf_depend_exit(1)
            else:
                kromsatel_args['reads_R2'].append(fpath)
            # end try

        elif opt in ('-u', '--reads-unpaired'):
            fpath = os.path.abspath(arg)
            try:
                _check_fastq_file(fpath)
            except ValueError as err:
                print_err(str(err))
                platf_depend_exit(1)
            else:
                kromsatel_args['reads_unpaired'].append(fpath)
            # end try

        elif opt in ('-d', '--db'):
            if not os.path.exists( '{}.nhr'.format(arg) ):
                print_err('Error: database `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                kromsatel_args['db_fpath'] = os.path.abspath(arg)
            # end if

        # Handle primers path
        elif opt in ('-p', '--primers'):
            if not os.path.exists(arg):
                print_err('Error: file `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            # end if
            kromsatel_args['primers_fpath'] = arg

        # Handle chunk size
        elif opt in ('-c', '--chunk-size'):
            try:
                kromsatel_args['chunk_size'] = int(arg)
                if kromsatel_args['chunk_size'] <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print_err('Error: chunk_size (option `{}`) must be positive integer number.'\
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
                print_err('Error: number of threads must be integer number > 0!')
                print_err(' And here is your value: `{}`'.format(arg))
                platf_depend_exit(1)
            # end try
            if kromsatel_args['n_thr'] > len(os.sched_getaffinity(0)):
                print_err('''\nWarning! You have specified {} threads to use
      although {} are available.'''.format(kromsatel_args['n_thr'], len(os.sched_getaffinity(0))))
                kromsatel_args['n_thr'] = len(os.sched_getaffinity(0))
                print_err('Switched to {} threads.\n'.format(kromsatel_args['n_thr']))
            # end if

        # Handle min major read length
        elif opt == '--am':
            try:
                kromsatel_args['min_len_major'] = int(arg)
                if kromsatel_args['min_len_major'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print_err('Invalid minimum major read length passed with option `{}`: {}'\
                    .format(opt, arg))
                print_err('It must be integer number > 0.')
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
                print_err('Invalid minimum minor read length passed with option `{}`: {}'\
                    .format(opt, arg))
                print_err('It must be integer number > 0.')
                platf_depend_exit(1)
            # end try

        elif opt in ('-o', '--outdir'):
            kromsatel_args['outdir'] = os.path.abspath(arg)
            if not os.path.isdir(kromsatel_args['outdir']):
                try:
                    os.makedirs(kromsatel_args['outdir'])
                except OSError as err:
                    print_err('Error: cannot create output directory.')
                    print_err(str(err))
                    platf_depend_exit(1)
                # end try
            # end if
        # end if
    # end for

    if not _mandatory_options_ok(kromsatel_args):
        platf_depend_exit(1)
    # end if

    if not _input_data_ok(kromsatel_args):
        platf_depend_exit(1)
    # end if

    if len(kromsatel_args['reads_R1']) != 0:
        kromsatel_args['paired_mode'] = True
        kromsatel_args['reads_R1'] = _get_first_element(kromsatel_args['reads_R1'])
        kromsatel_args['reads_R2'] = _get_first_element(kromsatel_args['reads_R2'])
        kromsatel_args['reads_unpaired'] = None
    elif len(kromsatel_args['reads_unpaired']) != 0:
        kromsatel_args['paired_mode'] = False
        kromsatel_args['reads_R1'] = None
        kromsatel_args['reads_R2'] = None
        kromsatel_args['reads_unpaired'] = _get_first_element(kromsatel_args['reads_unpaired'])
    # end if

    return kromsatel_args
# end def handle_cl_args


def _input_data_ok(kromsatel_args):

    mult_files_same_type = False
    mult_files_same_type = mult_files_same_type \
        or _mult_files_same_type(
            kromsatel_args['reads_R1'],
            file_type='forward (R1) reads'
        )
    mult_files_same_type = mult_files_same_type \
        or _mult_files_same_type(
            kromsatel_args['reads_R2'],
            file_type='reverse (R2) reads' 
        )
    mult_files_same_type = mult_files_same_type \
        or _mult_files_same_type(
            kromsatel_args['reads_unpaired'],
            file_type='unpaired reads'
        )
    if mult_files_same_type:
        return False
    # end if

    num_files_R1 = len(kromsatel_args['reads_R1'])
    num_files_R2 = len(kromsatel_args['reads_R2'])
    num_files_unpaired = len(kromsatel_args['reads_unpaired'])

    no_input_data = num_files_R1 == 0 \
                    and num_files_R2 == 0 \
                    and num_files_unpaired == 0
    if no_input_data:
        print_err('\nError: no input files have been specified.')
        return False
    # end if

    if num_files_R1 != num_files_R2:
        print_err('\nError: The number of files of forward reads ({} files) is not equal to \n\
    the number of files of reverse reads ({} files).'.format(num_files_R1, num_files_R2))
        return False
    # end if

    mixed_file_types = num_files_R1 != 0 \
                       and num_files_unpaired != 0
    if mixed_file_types:
        print_err('\nError: kromsatel cannot process "mixed" input data, \n\
    i.e. when both paired (R1/R2) and unpaired read files are specified.')
        print_err('You have specified {} pairs of paired files.'.format(num_files_R1))
        print_err('You have specified {} unpaired files.'.format(num_files_unpaired))
        return False
    # end if

    return True
# end def _input_data_ok


def _mandatory_options_ok(kromsatel_args):
    mandatory_options = ('-d/--db', '-p/--primers')
    mandatory_arg_names = ('db_fpath', 'primers_fpath')

    for opt, arg_name in zip(mandatory_options, mandatory_arg_names):
        if kromsatel_args[arg_name] is None:
            print_err('Error: option `{}` is mandatory.'.format(opt))
            print_err('Type {} -h to see help message.'.format(sys.argv[0]))
            return False
        # end if
    # end for

    return True
# end def _mandatory_options_ok


def _get_first_element(collection):
    return next(iter(collection))
# end def _get_first_element


def _mult_files_same_type(file_paths, file_type):
    if len(file_paths) > 1:
        print_err('\nError: multiple files of {} are not allowed!'.format(file_type))
        print_err('You have specified {} files of {}'.format(len(file_paths), file_type))
        for i, fpath in enumerate(file_paths):
            print_err('  {}. `{}`'.format(i+1, fpath))
        # end for
        return True
    # end if
    return False
# end def _mult_files_same_type


def _check_fastq_file(fpath):

    if not os.path.exists(fpath):
        raise ValueError('\nFile `{}` does not exist'.format(fpath))
    elif not _is_fastq(fpath):
        raise ValueError(
            """\nFile `{}` does not look like a fastq_file
    The script detects fastq files by the extention. It can process gzipped files.
    Permitted extentions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`""".format(fpath)
        )
    # end if
# end def _check_fastq_file


def _is_fastq(fpath):
    fastq_pattern = r'.*\.f(ast)?q(\.gz)?$'
    match_obj = re.match(fastq_pattern, fpath)
    return not match_obj is None
# end def _is_fastq
