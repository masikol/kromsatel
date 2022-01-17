
import os

import src.filesystem as fs


class Output:

    def __init__(self, outdir_path, output_prefix):
        self.outdir_path = outdir_path
        self.output_prefix = output_prefix
    # end def
# end class


class UnpairedOutput(Output):

    def __init__(self, outdir_path, output_prefix):

        super().__init__(outdir_path, output_prefix)

        self.major_outfpath = None
        self._set_major_outfpath()
        self.minor_outfpath = None
        self._set_minor_outfpath()
        self.uncertain_outfpath = None
        self._set_uncertain_outfpath()

        self._init_output()
    # end def __init__

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
    # end def init_output

    def _set_major_outfpath(self):
        suffix = 'major'
        self.major_outfpath = self._configure_outfpath(suffix)
    # end def _set_major_outfpath

    def _set_minor_outfpath(self):
        suffix = 'minor'
        self.minor_outfpath = self._configure_outfpath(suffix)
    # end def _set_minor_outfpath

    def _set_uncertain_outfpath(self):
        suffix = 'uncertain'
        self.uncertain_outfpath = self._configure_outfpath(suffix)
    # end def _set_uncertain_outfpath

    def _configure_outfpath(self, suffix):
        return os.path.join(
            self.outdir_path,
            '{}_{}.fastq.gz'.format(self.output_prefix, suffix)
        )
    # end def _configure_outfpath
# end class UnpairedOutput


class PairedOutput(Output):

    def __init__(self, outdir_path, output_prefix):
        super().__init__(outdir_path, output_prefix)
        self.sample_name = self._get_sample_name()

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
        self.major_frw_outfpath = self._configure_outfpath(suffix, forward=True)
        self.major_rvr_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def

    def _set_minor_outfpaths(self):
        suffix = 'minor'
        self.minor_frw_outfpath = self._configure_outfpath(suffix, forward=True)
        self.minor_rvr_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def

    def _set_uncertain_outfpaths(self):
        suffix = 'uncertain'
        self.uncertain_frw_outfpath = self._configure_outfpath(suffix, forward=True)
        self.uncertain_rvr_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def

    def _set_unpaired_outfpaths(self):
        suffix = 'unpaired'
        self.unpaired_frw_outfpath = self._configure_outfpath(suffix, forward=True)
        self.unpaired_rvr_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def

    def _get_sample_name(self):

        for direction in ('_R1_001', '_R2_001'):
            if direction in self.output_prefix:
                sample_name = self.output_prefix.replace(direction, '')
            # end if
        # end for
        return sample_name
    # end def

    def _configure_outfpath(self, suffix, forward=True):
        direction = 'R1_001' if forward else 'R2_001'
        return os.path.join(
            self.outdir_path,
            '{}_{}_{}.fastq.gz'.format(self.sample_name, direction, suffix)
        )
    # end def
# end class
