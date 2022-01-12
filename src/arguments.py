
import os

import src.blast
import src.filesystem as fs
from src.printing import print_err
from src.platform import platf_depend_exit


class KromsatelArgs:

    def __init__(self, argparse_args):
        self.argparse_args = argparse_args
        self._init_default_arguments()
        self._check_actual_arguments()
        self._set_actual_arguments()
    # end def

    def _init_default_arguments(self):
        # Input data
        self.forward_read_fpath = None
        self.reverse_read_fpath = None
        self.unpaired_read_fpath = None
        self.primers_fpath = None
        self.reference_fpath = None

        # Output
        self.outdir_path = os.path.join(
            os.getcwd(),
            'kromsatel_output'
        )
        self.min_len = 25 # bp

        # Computational resourses
        self.threads_num = 1 # thread

        # Advanced
        self.chunk_size = 1000 # reads
        self.blast_task = 'megablast'
        self.fixed_crop_len = 'auto'
        self.primer_ext_len = 5 # bp
        self.use_index = False

        # "Meta" parameter derived from arguments
        self.paired_mode = None
        self.tmp_dir_path = None
    # end def

    def __repr__(self):
        repr_str = """KromsatelArgs:
        paired_mode = {},
        forward_read_fpath = `{}`, reverse_read_fpath = `{}`,
        unpaired_read_fpath = `{}`,
        primers_fpath = `{}`,
        reference_fpath = `{}`,
        outdir_path = `{}`,
        min_len = {},
        threads_num = {},
        chunk_size = {},
        blast_task = {},
        fixed_crop_len = {},
        primer_ext_len = {},
        use_index = {}
        """.format(
            self.paired_mode,
            self.forward_read_fpath,
            self.reverse_read_fpath,
            self.unpaired_read_fpath,
            self.primers_fpath,
            self.reference_fpath,
            self.outdir_path,
            self.min_len,
            self.threads_num,
            self.chunk_size,
            self.blast_task,
            self.fixed_crop_len,
            self.primer_ext_len,
            self.use_index
        )
        return repr_str
    # end def

    def set_output(self, output):
        self.output = output
    # end def

    def set_database_path(self, df_fpath):
        self.db_fpath = df_fpath
    # end def

    def _check_actual_arguments(self):
        argument_checker = KromsatelArgumentChecker(self.argparse_args)
        argument_checker.check_arguments()
    # end def

    def _set_actual_arguments(self):
        self._set_reads_fpaths()
        self._set_primers_fpath()
        self._set_reference_fpath()
        self._set_outdpath()
        self._set_min_len()
        self._set_blast_task()
        self._set_fixed_crop_len()
        self._set_primer_ext_len()
        self._set_use_index()
    # end def

    def _set_reads_fpaths(self):
        self.paired_mode = not self.argparse_args.reads_R1 is None
        if self.paired_mode:
            self.forward_read_fpath = self.argparse_args.reads_R1
            self.reverse_read_fpath = self.argparse_args.reads_R2
        else:
            self.unpaired_read_fpath = self.argparse_args.reads_unpaired
        # end if
    # end def

    def _set_primers_fpath(self):
        self.primers_fpath = self.argparse_args.primers
    # end def

    def _set_reference_fpath(self):
        self.reference_fpath = self.argparse_args.reference
    # end def

    def _set_outdpath(self):
        if not self.argparse_args.outdir is None:
            self.outdir_path = self.argparse_args.outdir
        # end if
        self._create_tmp_directory()
    # end def

    def _create_tmp_directory(self):
        self.tmp_dir_path = os.path.join(
            self.outdir_path, 'tmp'
        )
        try:
            fs.create_dir(self.tmp_dir_path)
        except OSError as err:
            print_err('\nError: cannot create temporary directory `{}`' \
                .format(self.tmp_dir_path))
            print_err(err)
            print_err('The program will use the root of the output directory \n\
    as a storage of temporary files')
            self.tmp_dir_path = self.outdir_path
        # end try
    # end def

    def _set_min_len(self):
        if not self.argparse_args.min_len is None:
            min_len_string = self.argparse_args.min_len
            self.min_len = int(min_len_string)
        # end if
    # end def

    def _set_blast_task(self):
        if not self.argparse_args.blast_task is None:
            self.blast_task = self.argparse_args.blast_task
        # end if
    # end def

    def _set_fixed_crop_len(self):
        if not self.argparse_args.crop_len is None:
            self.fixed_crop_len = self.argparse_args.crop_len
        # end if
    # end def

    def _set_primer_ext_len(self):
        if not self.argparse_args.primer_5ext is None:
            self.primer_ext_len = self.argparse_args.primer_5ext
        # end if
    # end def

    def _set_use_index(self):
        if self.argparse_args.use_index is None:
            return
        # end if

        if self.argparse_args.use_index == 'auto':
            blast_task = self.argparse_args.blast_task
            if blast_task in src.blast.TASKS_SUPPORT_INDEXED_SEARCH:
                self.use_index = True
            else:
                self.use_index = False
            # end if
        else:
            if self.argparse_args.use_index == 'true':
                self.use_index = True
            else:
                self.use_index = False
            # end if
        # end if
    # end def
# end class


class KromsatelArgumentChecker:

    def __init__(self, argparse_args):
        self.argparse_args = argparse_args

    # end def

    def check_arguments(self):
        try:
            self._check_mandatory_args()

            self._check_reads_fpaths()
            self._check_primers_fpath()
            self._check_reference_fpath()
            self._check_outdpath()
            self._check_min_len()
            self._check_threads_num()
            self._check_chunk_size()
            self._check_blast_task()
            self._check_fixed_crop_len()
            self._check_primer_ext_len()
            self._check_use_index()
        except ArgumentError:
            platf_depend_exit(1)
        # end try
    # end def

    def _check_mandatory_args(self):
        mandatory_args = (
            self.argparse_args.primers,
            self.argparse_args.reference,
        )

        mandarory_args_descriptions = (
            '-p/--primers',
            '-r/--reference',
        )

        for arg, description in zip(mandatory_args, mandarory_args_descriptions):
            if arg is None:
                print_err('\nError: argument {} is mandatory'.format(description))
                raise ArgumentError
            # end if
        # end for
    # end def

    def _check_reads_fpaths(self):
        if not read_type_file_combnation_ok(self.argparse_args):
            raise ArgumentError
        # end if
        paired_mode = not self.argparse_args.reads_R1 is None
        if not self._reads_files_exist(paired_mode):
            raise ArgumentError
        # end if
    # end def

    def _reads_files_exist(self, paired_mode):

        file_paths_to_check_existance = None

        if paired_mode:
            file_paths_to_check_existance = (
                self.argparse_args.reads_R1,
                self.argparse_args.reads_R2,
            )
        else:
            file_paths_to_check_existance = (
                self.argparse_args.reads_unpaired,
            )
        # end if

        for file_path in file_paths_to_check_existance:
            if not os.path.exists(file_path):
                print('\nError: file `{}` does not exist'.format(file_path))
                return False
            # end if
        # end for
        return True
    # end def

    def _check_primers_fpath(self):
        if not os.path.exists(self.argparse_args.primers):
            print('\nError: file `{}` does not exist'.format(self.argparse_args.primers))
            raise ArgumentError
        # end if
    # end def

    def _check_reference_fpath(self):
        if not os.path.exists(self.argparse_args.reference):
            print('\nError: file `{}` does not exist'.format(self.argparse_args.reference))
            raise ArgumentError
        # end if
    # end def

    def _check_outdpath(self):
        if self.argparse_args.outdir is None:
            return
        # end if
        try:
            fs.create_dir(self.argparse_args.outdir)
        except OSError as err:
            print_err('\nError: cannot create directory `{}`'.format(self.argparse_args.outdir))
            print_err(str(err))
            raise ArgumentError
        # end if
    # end def

    def _check_min_len(self):
        if self.argparse_args.min_len is None:
            return
        # end if
        min_len_string = self.argparse_args.min_len
        try:
            check_int_string_gt0(min_len_string)
        except AtoiGreaterThanZeroError as err:
            print_err('\nError: invalid minimum length: `{}`'.format(min_len_string))
            print_err(str(err))
            raise ArgumentError
        # end try
    # end def

    def _check_threads_num(self):
        if self.argparse_args.threads is None:
            return
        # end if
        threads_num_string = self.argparse_args.threads
        try:
            check_int_string_gt0(threads_num_string)
        except AtoiGreaterThanZeroError as err:
            print_err('\nError: invalid number of threads: `{}`'.format(threads_num_string))
            print_err(str(err))
            raise ArgumentError
        # end try
    # end def

    def _check_chunk_size(self):
        if self.argparse_args.chunk_size is None:
            return
        # end if
        chunk_size_string = self.argparse_args.chunk_size
        try:
            check_int_string_gt0(chunk_size_string)
        except AtoiGreaterThanZeroError as err:
            print_err('\nError: invalid chunk size: `{}`'.format(chunk_size_string))
            print_err(str(err))
            raise ArgumentError
        # end try
    # end def

    def _check_blast_task(self):
        if self.argparse_args.blast_task is None:
            return
        # end if
        blast_task_argument = self.argparse_args.blast_task
        if not blast_task_argument in src.blast.BLAST_TASKS:
            print_err('\nError: invalid name of a blast task: `{}`.'.format(blast_task_argument))
            print_err('Allowed values: {}'.format(', '.join(src.blast.BLAST_TASKS)))
            raise ArgumentError
        # end if
    # end def

    def _check_fixed_crop_len(self):
        if self.argparse_args.crop_len is None:
            return
        # end if
        crop_len_string = self.argparse_args.crop_len
        auto_detect_crop_len = (crop_len_string == 'auto')
        if not auto_detect_crop_len:
            try:
                check_int_string_ge0(crop_len_string)
            except AtoiGreaterOrEqualToZeroError as err:
                print_err('\nError: invalid crop length: `{}`'.format(crop_len_string))
                print_err(str(err))
                print_err('Also, it may be `auto`.')
                raise ArgumentError
            # end try
        # end if
    # end def

    def _check_primer_ext_len(self):
        if self.argparse_args.primer_5ext is None:
            return
        # end if
        primer_ext_string = self.argparse_args.primer_5ext
        try:
            check_int_string_ge0(primer_ext_string)
        except AtoiGreaterOrEqualToZeroError as err:
            print_err('\nError: invalid size of primer coordinates extention: `{}`' \
                .format(primer_ext_string))
            print_err(str(err))
            raise ArgumentError
        # end try
    # end def

    def _check_use_index(self):
        if self.argparse_args.use_index is None:
            return
        # end if

        allowed_values = {'true', 'false', 'auto',}
        use_index_string = self.argparse_args.use_index

        if not use_index_string in allowed_values:
            print_err('\nError: invalid value of the `--use-index` option: `{}`' \
                .format(use_index_string))
            print_err('Allowed_values: {}'.format(', '.join(allowed_values)))
            raise ArgumentError
        # end if

        blast_task = self.argparse_args.blast_task

        blast_task_cant_use_index = use_index_string == 'true' \
                                    and not blast_task in src.blast.TASKS_SUPPORT_INDEXED_SEARCH
        if blast_task_cant_use_index:
            print_err('\nError: BLAST task {} cannot use indexed search'.format(blast_task))
            raise ArgumentError
        # end if
    # end def
# end class


class ArgumentError(Exception):
    pass
# end class


class AtoiGreaterThanZeroError(Exception):
    pass
# end class


class AtoiGreaterOrEqualToZeroError(Exception):
    pass
# end class


def check_int_string_gt0(string_value):
    try:
        int_value = int(string_value)
        if int_value < 1:
            raise ValueError
        # end if
    except ValueError:
        raise AtoiGreaterThanZeroError('This value must be integer > 0')
    # end try
# end def


def check_int_string_ge0(string_value):
    try:
        int_value = int(string_value)
        if int_value < 0:
            raise ValueError
        # end if
    except ValueError:
        raise AtoiGreaterThanZeroError('This value must be integer >= 0')
    # end try
# end def


def read_type_file_combnation_ok(argparse_args):

    read_pass_string = create_read_pass_string(argparse_args)

    no_input_data = (read_pass_string == 'fru')
    if no_input_data:
        print_err('\nError: no input data.')
        return False
    # end if

    not_enough_data_strings = (
        'Fru',
        'fRu',
    )
    if read_pass_string in not_enough_data_strings:
        print_err('\nError: not enough input data.')
        print_err('Files of forward reads passed:')
        enumerate_file_paths_to_stderr(argparse_args.reads_R1[0])
        print_err('---')
        print_err('Files of reverse reads passed:')
        enumerate_file_paths_to_stderr(argparse_args.reads_R2[0])
        print_err('---')
        return False
    # end if

    mixed_data_strings = (
        'fRU',
        'FrU',
        'FRU',
    )

    if read_pass_string in mixed_data_strings:
        print_err('\nError: kromsatel cannot process "mixed" input data, \n\
    i.e. when both paired (R1/R2) and unpaired read files are passed to it.')
        print_err('Files of forward reads passed:')
        enumerate_file_paths_to_stderr(argparse_args.reads_R1[0])
        print_err('---')
        print_err('Files of reverse reads passed:')
        enumerate_file_paths_to_stderr(argparse_args.reads_R2[0])
        print_err('---')
        print_err('Files of unpaired reads passed:')
        enumerate_file_paths_to_stderr(argparse_args.reads_unpaired[0])
        print_err('---')
        return False
    # end if

    return True
# end def


def enumerate_file_paths_to_stderr(file_paths):
    for i, fpath in enumerate(file_paths):
        print_err('  {}. `{}`'.format(i+1, fpath))
    # end for
# end def


def create_read_pass_string(argparse_args):

    forward_reads_passed  = not argparse_args.reads_R1       is None
    reverse_reads_passed  = not argparse_args.reads_R2       is None
    unpaired_reads_passed = not argparse_args.reads_unpaired is None

    forward_char  = 'F' if forward_reads_passed  else 'f'
    reverse_char  = 'R' if reverse_reads_passed  else 'r'
    unpaired_char = 'U' if unpaired_reads_passed else 'u'

    read_pass_string = '{}{}{}'.format(
        forward_char, reverse_char, unpaired_char
    )

    return read_pass_string
# end def
