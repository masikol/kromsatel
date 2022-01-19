
import os

import src.blast
import src.filesystem as fs
from src.printing import print_err
from src.fatal_errors import FatalError
from src.kromsatel_modes import KromsatelModes


class KromsatelArgs:

    def __init__(self, argparse_args):
        self.argparse_args = argparse_args
        self._init_default_arguments()
        self._check_actual_arguments()
        self._set_actual_arguments()
    # end def

    def _init_default_arguments(self):
        # Input data
        self.frw_read_fpath = None
        self.rvr_read_fpath = None
        self.long_read_fpath = None
        self.primers_fpath = None
        self.reference_fpath = None

        # Output
        self.outdir_path = os.path.join(
            os.getcwd(),
            'kromsatel_output'
        )
        self.split_output = False

        # Computational resourses
        self.threads_num = 1 # thread

        # Advanced
        self.min_len = 25 # bp
        self.chunk_size = 1000 # reads
        self.blast_task = 'megablast'
        self.fixed_crop_len = 'auto'
        self.primer_ext_len = 5 # bp
        self.use_index = False

        # "Meta" parameter derived from arguments
        self.kromsatel_mode = None
        self.tmp_dir_path = None
    # end def

    def __repr__(self):
        repr_str = 'KromsatelArgs:\n' \
        + 'kromsatel_mode = {}\n'   .format(self.kromsatel_mode) \
        + 'frw_read_fpath = `{}`\n' .format(self.frw_read_fpath) \
        + 'rvr_read_fpath = `{}`\n' .format(self.rvr_read_fpath) \
        + 'long_read_fpath = `{}`\n'.format(self.long_read_fpath) \
        + 'primers_fpath = `{}`\n'  .format(self.primers_fpath) \
        + 'reference_fpath = `{}`\n'.format(self.reference_fpath) \
        + 'outdir_path = `{}`\n'    .format(self.outdir_path) \
        + 'split_output = `{}`\n'   .format(self.split_output) \
        + 'min_len = {}\n'          .format(self.min_len) \
        + 'threads_num = {}\n'      .format(self.threads_num) \
        + 'chunk_size = {}\n'       .format(self.chunk_size) \
        + 'blast_task = {}\n'       .format(self.blast_task) \
        + 'fixed_crop_len = {}\n'   .format(self.fixed_crop_len) \
        + 'primer_ext_len = {}\n'   .format(self.primer_ext_len) \
        + 'use_index = {}\n'        .format(self.use_index)
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
        self._set_split_output()
        self._set_min_len()
        self._set_threads_num()
        self._set_blast_task()
        self._set_fixed_crop_len()
        self._set_primer_ext_len()
        self._set_use_index()
    # end def

    def _set_reads_fpaths(self):
        self.kromsatel_mode = _detect_kromsatel_mode(self.argparse_args)
        if self.kromsatel_mode == KromsatelModes.IlluminaPE:
            self.frw_read_fpath = self.argparse_args.reads_R1
            self.rvr_read_fpath = self.argparse_args.reads_R2
        elif self.kromsatel_mode == KromsatelModes.Nanopore:
            self.long_read_fpath = self.argparse_args.reads_long
        elif self.kromsatel_mode == KromsatelModes.IlluminaSE:
            self.frw_read_fpath = self.argparse_args.reads_R1
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
        except FatalError as err:
            print_err('\nError: cannot create temporary directory `{}`' \
                .format(self.tmp_dir_path))
            print_err(err)
            print_err('The program will use the root of the output directory' \
                ' as a storage of temporary files')
            self.tmp_dir_path = self.outdir_path
        # end try
    # end def

    def _set_split_output(self):
         self.split_output = self.argparse_args.split_output
    # end def

    def _set_min_len(self):
        if not self.argparse_args.min_len is None:
            min_len_string = self.argparse_args.min_len
            self.min_len = int(min_len_string)
        # end if
    # end def

    def _set_threads_num(self):
        if not self.argparse_args.threads is None:
            num_threads_string = self.argparse_args.threads
            self.threads_num = int(num_threads_string)
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
                error_msg = '\nError: argument {} is mandatory'.format(description)
                raise FatalError(error_msg)
            # end if
        # end for
    # end def

    def _check_reads_fpaths(self):
        try:
            _check_file_type_combination(self.argparse_args)
        except _InvalidFileCombinationError as err:
            raise FatalError(str(err))
        # end try

        kromsatel_mode = _detect_kromsatel_mode(self.argparse_args)
        try:
            self._reads_files_exist(kromsatel_mode)
        except FileNotFoundError as err:
            raise FatalError(str(err))
        # end try
    # end def

    def _reads_files_exist(self, kromsatel_mode):

        file_paths_to_check_existance = None

        if kromsatel_mode == KromsatelModes.IlluminaPE:
            file_paths_to_check_existance = (
                self.argparse_args.reads_R1,
                self.argparse_args.reads_R2,
            )
        elif kromsatel_mode == KromsatelModes.Nanopore:
            file_paths_to_check_existance = (
                self.argparse_args.reads_long,
            )
        elif kromsatel_mode == KromsatelModes.IlluminaSE:
            file_paths_to_check_existance = (
                self.argparse_args.reads_R1,
            )
        # end if

        non_extant_file_paths = list()

        for file_path in file_paths_to_check_existance:
            if not os.path.exists(file_path):
                non_extant_file_paths.append(file_path)
            # end if
        # end for

        if len(non_extant_file_paths) != 0:
            error_msg = '\nError: following files do not exits:\n'
            for i, file_path in enumerate(non_extant_file_paths):
                error_msg += '  {}. `{}`\n'.format(i+1, file_path)
            # end for
            raise FileNotFoundError(error_msg)
        # end if
    # end def

    def _check_primers_fpath(self):
        if not os.path.exists(self.argparse_args.primers):
            error_msg = '\nError: file `{}` does not exist' \
                .format(self.argparse_args.primers)
            raise FatalError(error_msg)
        # end if
    # end def

    def _check_reference_fpath(self):
        if not os.path.exists(self.argparse_args.reference):
            error_msg = '\nError: file `{}` does not exist' \
                .format(self.argparse_args.reference)
            raise FatalError(error_msg)
        # end if
    # end def

    def _check_outdpath(self):
        if self.argparse_args.outdir is None:
            return
        # end if
        try:
            fs.create_dir(self.argparse_args.outdir)
        except FatalErrors as err:
            error_msg = '\nError: cannot create directory `{}`:\n {}' \
                .format(self.argparse_args.outdir, err)
            raise FatalError(error_msg)
        # end if
    # end def

    def _check_min_len(self):
        if self.argparse_args.min_len is None:
            return
        # end if
        min_len_string = self.argparse_args.min_len
        try:
            _check_int_string_gt0(min_len_string)
        except _AtoiGreaterThanZeroError as err:
            error_msg = '\nError: invalid minimum length: `{}`:\n {}' \
                .format(min_len_string, err)
            raise FatalError(error_msg)
        # end try
    # end def

    def _check_threads_num(self):
        if self.argparse_args.threads is None:
            return
        # end if
        threads_num_string = self.argparse_args.threads
        try:
            _check_int_string_gt0(threads_num_string)
        except _AtoiGreaterThanZeroError as err:
            error_msg = '\nError: invalid number of threads: `{}`\n {}' \
                .format(threads_num_string, err)
            raise FatalError(error_msg)
        # end try
    # end def

    def _check_chunk_size(self):
        if self.argparse_args.chunk_size is None:
            return
        # end if
        chunk_size_string = self.argparse_args.chunk_size
        try:
            _check_int_string_gt0(chunk_size_string)
        except _AtoiGreaterThanZeroError as err:
            error_msg = '\nError: invalid chunk size: `{}`\n {}' \
                .format(chunk_size_string, err)
            raise FatalError(error_msg)
        # end try
    # end def

    def _check_blast_task(self):
        if self.argparse_args.blast_task is None:
            return
        # end if
        blast_task_argument = self.argparse_args.blast_task
        if not blast_task_argument in src.blast.BLAST_TASKS:
            error_msg = '\nError: invalid name of a blast task: `{}`.' \
                'Allowed values: {}' \
                .format(blast_task_argument, ', '.join(src.blast.BLAST_TASKS))
            raise FatalError(error_msg)
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
                _check_int_string_ge0(crop_len_string)
            except _AtoiGreaterOrEqualToZeroError as err:
                error_msg = '\nError: invalid crop length: `{}`\n  {}' \
                    '  Also, it may be `auto`.' \
                    .format(crop_len_string, err)
                raise FatalError(error_msg)
            # end try
        # end if
    # end def

    def _check_primer_ext_len(self):
        if self.argparse_args.primer_5ext is None:
            return
        # end if
        primer_ext_string = self.argparse_args.primer_5ext
        try:
            _check_int_string_ge0(primer_ext_string)
        except _AtoiGreaterOrEqualToZeroError as err:
            error_msg = '\nError: invalid size of primer coordinates extention: `{}`\n  {}' \
                .format(primer_ext_string, err)
            raise FatalError(error_msg)
        # end try
    # end def

    def _check_use_index(self):
        if self.argparse_args.use_index is None:
            return
        # end if

        allowed_values = {'true', 'false', 'auto',}
        use_index_string = self.argparse_args.use_index

        if not use_index_string in allowed_values:
            error_msg = '\nError: invalid value of the `--use-index` option: `{}`' \
                        'Allowed_values: {}' \
                            .format(use_index_string, ', '.join(allowed_values))
            raise FatalError(error_msg)
        # end if

        blast_task = self.argparse_args.blast_task

        blast_task_cant_use_index = use_index_string == 'true' \
                                    and not blast_task in src.blast.TASKS_SUPPORT_INDEXED_SEARCH
        if blast_task_cant_use_index:
            error_msg = '\nError: BLAST task {} cannot use indexed search' \
                .format(blast_task)
            raise FatalError(error_msg)
        # end if
    # end def
# end class


class _AtoiGreaterThanZeroError(Exception):
    pass
# end class


class _AtoiGreaterOrEqualToZeroError(Exception):
    pass
# end class


class _InvalidFileCombinationError(Exception):
    pass
# end class


def _detect_kromsatel_mode(argparse_args):
    read_pass_string = _create_read_pass_string(argparse_args)

    if read_pass_string == 'FRl':
        return KromsatelModes.IlluminaPE
    elif read_pass_string == 'frL':
        return KromsatelModes.Nanopore
    elif read_pass_string == 'Frl':
        return KromsatelModes.IlluminaSE
    # end if

    # Execution should not reach here
    error_msg = '\nInternal error. Please, contact the developer' \
        ' and tell him abouth this error.\n' \
        'Error description: "kromsatel mode error in _detect_kromsatel_mode".'
    raise FatalError(error_msg)
# end def


def _check_int_string_gt0(string_value):
    try:
        int_value = int(string_value)
        if int_value < 1:
            raise ValueError
        # end if
    except ValueError:
        raise _AtoiGreaterThanZeroError('This value must be integer > 0')
    # end try
# end def


def _check_int_string_ge0(string_value):
    try:
        int_value = int(string_value)
        if int_value < 0:
            raise ValueError
        # end if
    except ValueError:
        raise _AtoiGreaterThanZeroError('This value must be integer >= 0')
    # end try
# end def


def _check_file_type_combination(argparse_args):

    read_pass_string = _create_read_pass_string(argparse_args)

    no_input_data = (read_pass_string == 'frl')
    if no_input_data:
        raise _InvalidFileCombinationError('\nError: no input data.')
    # end if

    mixed_data_strings = (
        'fRL',
        'FrL',
        'FRL',
    )

    if read_pass_string in mixed_data_strings:
        error_msg = _make_mixed_data_error_msg(
            argparse_args.reads_R1,
            argparse_args.reads_R2,
            argparse_args.reads_long
        )
        raise _InvalidFileCombinationError(error_msg)
    # end if
# end def


def _make_mixed_data_error_msg(frw_fpath, rvr_fpath, long_fpath):
    error_msg = '\nError: the program cannot process "mixed" input data,\n' \
                '  i.e. when both paired (R1/R2) and long read files are passed to it.\n' \
                'File of forward reads passed:\n  `{}`\n' \
                'File of reverse reads passed:\n  `{}`\n' \
                'File of long reads passed:\n  `{}`\n' \
                    .format(frw_fpath, rvr_fpath, long_fpath)
    return error_msg
# end def


def _create_read_pass_string(argparse_args):

    frw_reads_passed  = not argparse_args.reads_R1   is None
    rvr_reads_passed  = not argparse_args.reads_R2   is None
    long_reads_passed = not argparse_args.reads_long is None

    frw_char  = 'F' if frw_reads_passed  else 'f'
    rvr_char  = 'R' if rvr_reads_passed  else 'r'
    long_char = 'L' if long_reads_passed else 'l'

    read_pass_string = '{}{}{}'.format(
        frw_char, rvr_char, long_char
    )

    return read_pass_string
# end def
