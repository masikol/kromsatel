
import os

import src.filesystem as fs


class Output:

    def __init__(self, outdir_path, output_prefix):
        self.outdir_path = outdir_path
        self.output_prefix = output_prefix
    # end def
# end class


class SimpleUnpairedOutput(Output):

    def __init__(self, outdir_path, output_prefix):

        super().__init__(outdir_path, output_prefix)

        self.outfpath = None
        self._set_outfpath()
        self._init_output()
    # end def

    def _set_outfpath(self):
        suffix = 'cleaned'
        self.outfpath = _configure_unpaired_outfpath(
            self.outdir_path,
            self.output_prefix,
            suffix
        )
    # end def

    def _init_output(self):
        fs.create_dir(self.outdir_path)
        fs.init_file(self.outfpath)
    # end def
# end class


class SimplePairedOutput(Output):

    def __init__(self, outdir_path, output_prefix):

        super().__init__(outdir_path, output_prefix)
        self.sample_name = _get_sample_name(self.output_prefix)

        self.frw_outfpath = None
        self.rvr_outfpath = None
        self._set_normal_outfpaths()
        self.frw_unpaired_outfpath = None
        self.rvr_unpaired_outfpath = None
        self._set_unpaired_outfpaths()

        self._init_output()
    # end def

    def _set_normal_outfpaths(self):
        suffix = 'cleaned'
        self.frw_outfpath = _configure_paired_outfpath(
            self.outdir_path,
            self.sample_name,
            suffix,
            forward=True
        )
        self.rvr_outfpath = _configure_paired_outfpath(
            self.outdir_path,
            self.sample_name,
            suffix,
            forward=False
        )
    # end def

    def _set_unpaired_outfpaths(self):
        suffix = 'cleaned_unpaired'
        self.frw_unpaired_outfpath = _configure_paired_outfpath(
            self.outdir_path,
            self.sample_name,
            suffix,
            forward=True
        )
        self.rvr_unpaired_outfpath = _configure_paired_outfpath(
            self.outdir_path,
            self.sample_name,
            suffix,
            forward=False
        )
    # end def

    def _init_output(self):

        fs.create_dir(self.outdir_path)

        output_fpaths = (
            self.frw_outfpath,
            self.rvr_outfpath,
            self.frw_unpaired_outfpath,
            self.rvr_unpaired_outfpath,
        )
        for outfpath in output_fpaths:
            fs.init_file(outfpath)
        # end for
    # end def
# end class


class SplitUnpairedOutput(Output):

    def __init__(self, outdir_path, output_prefix):

        super().__init__(outdir_path, output_prefix)

        self.major_outfpath = None
        self._set_major_outfpath()
        self.minor_outfpath = None
        self._set_minor_outfpath()
        self.uncertain_outfpath = None
        self._set_uncertain_outfpath()

        self._init_output()
    # end def

    def _init_output(self):

        fs.create_dir(self.outdir_path)

        output_fpaths = (
            self.major_outfpath,
            self.minor_outfpath,
            self.uncertain_outfpath,
        )
        for outfpath in output_fpaths:
            fs.init_file(outfpath)
        # end for
    # end def

    def _set_major_outfpath(self):
        suffix = 'major'
        self.major_outfpath = \
            _configure_unpaired_outfpath(
                self.outdir_path,
                self.output_prefix,
                suffix
            )
    # end def

    def _set_minor_outfpath(self):
        suffix = 'minor'
        self.minor_outfpath = \
            _configure_unpaired_outfpath(
                self.outdir_path,
                self.output_prefix,
                suffix
            )
    # end def

    def _set_uncertain_outfpath(self):
        suffix = 'uncertain'
        self.uncertain_outfpath = \
            _configure_unpaired_outfpath(
                self.outdir_path,
                self.output_prefix,
                suffix
            )
    # end def
# end class


class SplitPairedOutput(Output):

    def __init__(self, outdir_path, output_prefix):
        super().__init__(outdir_path, output_prefix)
        self.sample_name = _get_sample_name(self.output_prefix)

        self.major_frw_outfpath = None
        self.major_rvr_outfpath = None
        self._set_major_outfpaths()
        self.minor_frw_outfpath = None
        self.minor_rvr_outfpath = None
        self._set_minor_outfpaths()
        self.uncertain_frw_outfpath = None
        self.uncertain_rvr_outfpath = None
        self._set_uncertain_outfpaths()
        self.unpaired_frw_outfpath = None
        self.unpaired_rvr_outfpath = None
        self._set_unpaired_outfpaths()

        self._init_output()
    # end def

    def _init_output(self):
        fs.create_dir(self.outdir_path)

        output_fpaths = (
            self.major_frw_outfpath,
            self.major_rvr_outfpath,
            self.minor_frw_outfpath,
            self.minor_rvr_outfpath,
            self.uncertain_frw_outfpath,
            self.uncertain_rvr_outfpath,
            self.unpaired_frw_outfpath,
            self.unpaired_rvr_outfpath,
        )
        for outfpath in output_fpaths:
            fs.init_file(outfpath)
        # end for
    # end def

    def _set_major_outfpaths(self):
        suffix = 'major'
        self.major_frw_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=True
            )
        self.major_rvr_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=False
            )
    # end def

    def _set_minor_outfpaths(self):
        suffix = 'minor'
        self.minor_frw_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=True
            )
        self.minor_rvr_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=False
            )
    # end def

    def _set_uncertain_outfpaths(self):
        suffix = 'uncertain'
        self.uncertain_frw_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=True
            )
        self.uncertain_rvr_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=False
            )
    # end def

    def _set_unpaired_outfpaths(self):
        suffix = 'unpaired'
        self.unpaired_frw_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=True
            )
        self.unpaired_rvr_outfpath = \
            _configure_paired_outfpath(
                self.outdir_path,
                self.sample_name,
                suffix,
                forward=False
            )
    # end def
# end class


def _get_sample_name(output_prefix):

    for direction in ('_R1_001', '_R2_001'):
        if direction in output_prefix:
            sample_name = output_prefix.replace(direction, '')
        # end if
    # end for
    return sample_name
# end def


def _configure_unpaired_outfpath(outdir_path, output_prefix, suffix):
    return os.path.join(
        outdir_path,
        '{}_{}.fastq.gz'.format(output_prefix, suffix)
    )
# end def


def _configure_paired_outfpath(outdir_path, sample_name, suffix, forward=True):
    direction = 'R1_001' if forward else 'R2_001'
    return os.path.join(
        outdir_path,
        '{}_{}_{}.fastq.gz'.format(sample_name, direction, suffix)
    )
# end def
