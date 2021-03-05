# -*- coding: utf-8 -*-

import os
import subprocess as sp

from src.platform import platf_depend_exit


def check_program(program):
    # Check if program is in PATH
    pathdirs = os.environ['PATH'].split(os.pathsep)
    utility_found = False
    for directory in pathdirs:
        if os.path.exists(directory) and program in os.listdir(directory):
            utility_found = True
            break
        # end if
    # end for
    if not utility_found:
        print('  Attention!\n`{}` from BLAST+ toolkit is not installed.'.format(program))
        print("""If this error still occures although you have installed everything
 -- make sure that this program is added to PATH)""")
        platf_depend_exit(1)
    # end if
# end def check_blastn


def check_blast_db(db_fpath):
    # Function checks if database passed to kromsatel is usable.
    # It checkis it with blastdbcmd.

    program = 'blastdbcmd'
    check_program(program)

    check_cmd = '{} -info -db {}'.format(program, db_fpath)

    # Check database with blastdbcmd
    pipe = sp.Popen(check_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError: blast database at `{}` is not usable.'.format(db_fpath))
        print('Please, make it again.')
        print(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(pipe.returncode)
    # end if
# end def check_blast_db


def configure_blastn_cmd(query_fpath, db_fpath):
    # Function creates command for launching blastn.
    # This command will be always the same, so let's create it once.
    #
    # :param query_fpath: path to query fasta file;
    # :type query_fpath: str;
    # :param db_fpath: path to database;
    # :type db_fpath: str;
    #
    # Returns command which will launch aligning with blastn (str).

    outfmt = '6 qseqid sseqid sstrand length qlen slen qstart qend sstart send'

    # Configure command line
    blast_cmd = "blastn -query {} \
    -db {} \
    -task {} \
    -evalue 1e-3 \
    -outfmt '{}'".format(query_fpath,
                         db_fpath,
                         'dc-megablast',
                         outfmt
    )

    return blast_cmd
# end def configure_blastn_cmd


def blast_align(blast_cmd):
    # Function alignes reads against fragments using Discontiguous Megablast.
    #
    # :param blast_cmd: command to execute;
    # :type blast_cmd: str;
    #
    # Returns tabular string if format "outfmt 6" returned by blastn.

    # Launch blastn
    pipe = sp.Popen(blast_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError while aligning a sequence against amplion database')
        print(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(pipe.returncode)
    # end if

    return stdout_stderr[0].decode('utf-8')
# end def blast_align
