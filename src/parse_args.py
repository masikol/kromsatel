
import re
import os
import sys
import getopt

import src.blast
import src.filesystem as fs
from src.printing import print_err
from src.platform import platf_depend_exit


def handle_cl_args():
    # Function handles command line arguments.

    # Get arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hv1:2:u:p:r:t:c:m:k:o:',
            [
                'help', 'version',
                'reads-R1=', 'reads-R2=', 'reads-unpaired=',
                'primers=', 'reference=',
                'threads=', 'chunk-size=', 'min-len=', 'blast-task=',
                'crop-len='
                'outdir=',
            ]
        )
    except getopt.GetoptError as opt_err:
        print_err( str(opt_err) )
        print_err('See help (`-h` option)')
        platf_depend_exit(2)
    # end try

    kromsatel_args = {
        'reads_R1': list(),
        'reads_R2': list(),
        'reads_unpaired': list(),
        'primers_fpath': None,     # path to file of primers to remove
        'reference_fpath': None,   # path to database
        'n_thr': 1,                # number of threads
        'chunk_size': 1000,        # reads
        'min_len': 25,             # minimum length of an output read (bp)
        'blast_task': 'megablast', # variant of BLAST algorithm
        'fixed_crop_len': 27,      # size of sequence to crop from reads of non-specific amplicons
        'outdir': os.path.join(    # output directory
            os.getcwd(),
            'kromsatel_output'
        ),
        'paired_mode': False,      # True if PE, False if unpaired
        'use_index': True,         # True if blast should perform indexed search
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

        # Handle primers path
        elif opt in ('-p', '--primers'):
            if not os.path.exists(arg):
                print_err('\nError: file `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            # end if
            kromsatel_args['primers_fpath'] = arg

        elif opt in ('-r', '--reference'):
            if not os.path.exists('{}'.format(arg)):
                print_err('\nError: file `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                kromsatel_args['reference_fpath'] = os.path.abspath(arg)
            # end if

        # Handle chunk size
        elif opt in ('-c', '--chunk-size'):
            try:
                kromsatel_args['chunk_size'] = int(arg)
                if kromsatel_args['chunk_size'] <= 0:
                    raise ValueError
                # end if
            except ValueError:
                print_err('\nError: chunk_size (option `{}`) must be positive integer number.'\
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
                print_err('\nError: number of threads must be integer number > 0!')
                print_err(' And here is your value: `{}`'.format(arg))
                platf_depend_exit(1)
            # end try
            if kromsatel_args['n_thr'] > len(os.sched_getaffinity(0)):
                print_err('''\nWarning! You have specified {} threads to use
      although {} are available.'''.format(kromsatel_args['n_thr'], len(os.sched_getaffinity(0))))
                kromsatel_args['n_thr'] = len(os.sched_getaffinity(0))
                print_err('Switched to {} threads.\n'.format(kromsatel_args['n_thr']))
            # end if

        # Handle min minor read length
        elif opt in ('-m', '--min-len'):
            try:
                kromsatel_args['min_len'] = int(arg)
                if kromsatel_args['min_len'] < 1:
                    raise ValueError
                # end if
            except ValueError:
                print_err('\nInvalid value passed with the option `{}`: {}'\
                    .format(opt, arg))
                print_err('It must be integer number > 0.')
                platf_depend_exit(1)
            # end try

        elif opt in ('-k', '--blast-task'):
            if not arg in src.blastn.BLAST_TASKS:
                print_err('\nError: invalid name of a blast task: {}.'.format(arg))
                print_err('Allowed values: {}'.format(', '.join(src.blastn.BLAST_TASKS)))
                platf_depend_exit(1)
            # end if
            kromsatel_args['blast_task'] = arg

        elif opt == '--crop-len':
            try:
                kromsatel_args['fixed_crop_len'] = int(arg)
                if kromsatel_args['fixed_crop_len'] < 0:
                    raise ValueError
                # end if
            except ValueError:
                print_err('\nInvalid value passed with the option `{}`: {}'\
                    .format(opt, arg))
                print_err('It must be integer number >= 0.')
                platf_depend_exit(1)
            # end try

        elif opt in ('-o', '--outdir'):
            kromsatel_args['outdir'] = os.path.abspath(arg)
            if not os.path.isdir(kromsatel_args['outdir']):
                try:
                    os.makedirs(kromsatel_args['outdir'])
                except OSError as err:
                    print_err('\nError: cannot create output directory.')
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

    kromsatel_args['use_index'] = src.blast.if_use_index(kromsatel_args['blast_task'])

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

    kromsatel_args['tmp_dir'] = os.path.join(kromsatel_args['outdir'], 'tmp')
    fs.create_dir(kromsatel_args['tmp_dir'])

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
    mandatory_options = ('-r/--reference', '-p/--primers')
    mandatory_arg_names = ('reference_fpath', 'primers_fpath')

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
