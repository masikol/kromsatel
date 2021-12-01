# -*- coding: utf-8 -*-

import re
import os
import gzip

from src.platform import platf_depend_exit

# For opening plain text and gzipped files
OPEN_FUNCS = (open, gzip.open)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode('utf-8').strip()  # format gzipped line
)


def make_outfpath(fq_fpath, outdir):
    # Function configures path to output file.
    #
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath: str;
    # :param outdir: path to output directory;
    # :type outdir: str;
    #
    # Returns path to output file (str).

    bname = os.path.basename(fq_fpath)

    try:
        extention = re.search(r'(\.f(ast)?q(\.gz)?)', bname).group(1)
    except AttributeError as err:
        print( str(err) )
        print('Error -3. Please, contact the developer.')
        platf_depend_exit(-3)
    # end try

    name_with_no_ext = bname.replace(extention, '')

    return os.path.join(
        outdir,
        '{}_cleaned.fastq'.format(name_with_no_ext)
    )
# end def make_outfpath


def rm_query_file(query_fpath):
    # Function removes temporary query file.
    #
    # :param query_fpath: path to query file;
    # :type query_fpath: str;

    try:
        os.unlink(query_fpath)
    except OSError as oserr:
        print('Warning: Cannot remove temporary file `{}`'.format(query_fpath))
        print( str(oserr) )
    # end try
# end rm_query_file
