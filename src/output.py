
import os

import src.filesystem as fs


class UnpairedOutput:

    def __init__(self, args):
        self.outdir = args['outdir']
        self.input_basename = fs.rm_fastq_extention(
            os.path.basename(args['reads_unpaired'])
        )

        self.major_outfpath = None
        self._set_major_outfpath()
        self.minor_outfpath = None
        self._set_minor_outfpath()
        self.trash_outfpath = None
        self._set_trash_outfpath()

        self._init_output()
    # end def __init__

    def _init_output(self):

        fs.create_dir(self.outdir)

        output_fpaths = (
            self.major_outfpath,
            self.minor_outfpath,
            self.trash_outfpath,
        )
        for outfpath in output_fpaths:
            fs.init_file(outfpath)
        # end for
    # end def init_output

    def _set_major_outfpath(self):
        suffix = 'cleaned'
        self.major_outfpath = self._configure_outfpath(suffix)
    # end def _set_major_outfpath

    def _set_minor_outfpath(self):
        suffix = 'minor_fragments'
        self.minor_outfpath = self._configure_outfpath(suffix)
    # end def _set_minor_outfpath

    def _set_trash_outfpath(self):
        suffix = 'trash'
        self.trash_outfpath = self._configure_outfpath(suffix)
    # end def _set_trash_outfpath

    def _configure_outfpath(self, suffix):
        return os.path.join(
            self.outdir,
            '{}_{}.fastq'.format(self.input_basename, suffix)
        )
    # end def _configure_outfpath
# end class UnpairedOutput


class PairedOutput:

    def __init__(self, args):
        self.outdir = args['outdir']
        self.input_basename = fs.rm_fastq_extention(
            os.path.basename(args['reads_R1'])
        )
        self.sample_name = self._get_sample_name()

        self.forward_major_outfpath = None
        self.reverse_major_outfpath = None
        self._set_major_outfpaths()
        self.forward_minor_outfpath = None
        self.reverse_minor_outfpath = None
        self._set_minor_outfpaths()
        self.forward_trash_outfpath = None
        self.reverse_trash_outfpath = None
        self._set_trash_outfpaths()

        self._init_output()
    # end def __init__

    def _init_output(self):
        fs.create_dir(self.outdir)

        output_fpaths = (
            self.forward_major_outfpath,
            self.reverse_major_outfpath,
            self.forward_minor_outfpath,
            self.reverse_minor_outfpath,
            self.forward_trash_outfpath,
            self.reverse_trash_outfpath,
        )
        for outfpath in output_fpaths:
            fs.init_file(outfpath)
        # end for
    # end def init_output

    def _set_major_outfpaths(self):
        suffix = 'cleaned'
        self.forward_major_outfpath = self._configure_outfpath(suffix, forward=True)
        self.reverse_major_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def _set_major_outfpaths

    def _set_minor_outfpaths(self):
        suffix = 'minor_fragments'
        self.forward_minor_outfpath = self._configure_outfpath(suffix, forward=True)
        self.reverse_minor_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def _set_minor_outfpaths

    def _set_trash_outfpaths(self):
        suffix = 'trash'
        self.forward_trash_outfpath = self._configure_outfpath(suffix, forward=True)
        self.reverse_trash_outfpath = self._configure_outfpath(suffix, forward=False)
    # end def _set_trash_outfpaths

    def _get_sample_name(self):

        for direction in ('_R1_001', '_R2_001'):
            if direction in self.input_basename:
                sample_name = self.input_basename.replace(direction, '')
            # end if
        # end for
        return sample_name
    # end def _get_sample_from_fastq_basename

    def _configure_outfpath(self, suffix, forward=True):
        direction = 'R1_001' if forward else 'R2_001'
        return os.path.join(
            self.outdir,
            '{}_{}_{}.fastq'.format(self.sample_name, direction, suffix)
        )
    # end def _configure_outfpath
# end class UnpairedOutput
