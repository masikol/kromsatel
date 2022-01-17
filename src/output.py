
import os

import src.filesystem as fs


class UnpairedOutput:

    def __init__(self, kromsatel_args):
        self.outdir = kromsatel_args.outdir_path
        self.input_basename = fs.rm_fastq_extention(
            os.path.basename(kromsatel_args.unpaired_read_fpath)
        )

        self.major_outfpath = None
        self._set_major_outfpath()
        self.minor_outfpath = None
        self._set_minor_outfpath()
        self.uncertain_outfpath = None
        self._set_uncertain_outfpath()

        self._init_output()
    # end def __init__

    def _init_output(self):

        fs.create_dir(self.outdir)

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
            self.outdir,
            '{}_{}.fastq.gz'.format(self.input_basename, suffix)
        )
    # end def _configure_outfpath
# end class UnpairedOutput


class PairedOutput:

    def __init__(self, kromsatel_args):
        self.outdir = kromsatel_args.outdir_path
        self.input_basename = fs.rm_fastq_extention(
            os.path.basename(kromsatel_args.frw_read_fpath)
        )
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
        fs.create_dir(self.outdir)

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
            if direction in self.input_basename:
                sample_name = self.input_basename.replace(direction, '')
            # end if
        # end for
        return sample_name
    # end def

    def _configure_outfpath(self, suffix, forward=True):
        direction = 'R1_001' if forward else 'R2_001'
        return os.path.join(
            self.outdir,
            '{}_{}_{}.fastq.gz'.format(self.sample_name, direction, suffix)
        )
    # end def
# end class
