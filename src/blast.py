
import os
import json
import subprocess as sp

import src.fastq
import src.filesystem as fs
from src.printing import print_err
from src.platform import platf_depend_exit


BLAST_TASKS = (
    'megablast',
    'dc-megablast',
    'blastn',
)


def get_blastplus_dependencies(krosatel_args):

    dependencies = [
        'blastn',
        'makeblastdb'
    ]

    if krosatel_args['use_index']:
        dependencies.append('makembindex')
    # end if

    return dependencies
# end def


def if_use_index(blast_task):

    # `megablast` and `blastn` support indexed searches
    # `dc-megablast` does not
    megablast, discomegablast, blastn = range(3)
    blast_tasks_use_index = (
        BLAST_TASKS[megablast],
        BLAST_TASKS[blastn],
    )
    blast_tasks_dont_use_index = (
        BLAST_TASKS[discomegablast],
    )

    if blast_task in blast_tasks_use_index:
        return True
    elif blast_task in blast_tasks_dont_use_index:
        return False
    else:
        raise ValueError(
            'Invalid name of a blast task: {}.\n  Allowed values: {}'.format(
                blast_task,
                ', '.join(BLAST_TASKS)
            )
        )
    # end if
# end def


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


# def check_blast_db(db_fpath):
#     # Function checks if database passed to kromsatel is usable.
#     # It checkis it with blastdbcmd.

#     program = 'blastdbcmd'
#     check_program(program)

#     check_cmd = '{} -info -db {}'.format(program, db_fpath)

#     # Check database with blastdbcmd
#     pipe = sp.Popen(check_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
#     stdout_stderr = pipe.communicate()

#     if pipe.returncode != 0:
#         print('\nError: blast database at `{}` is not usable.'.format(db_fpath))
#         print('Please, make it again.')
#         print(stdout_stderr[1].decode('utf-8'))
#         platf_depend_exit(pipe.returncode)
#     # end if
# # end def check_blast_db


def _configure_makeblastdb_cmd(fasta_fpath, db_fpath):
    makeblastdb_cmd = ' '.join(
        [
            'makeblastdb',
            '-in {} -input_type fasta'.format(fasta_fpath),
            '-dbtype nucl',
            '-parse_seqids',
            '-out {}'.format(db_fpath),
        ]
    )

    return makeblastdb_cmd
# end def


def _configure_makembindex_cmd(db_fpath):
    makembindex_cmd = ' '.join(
        [
            'makembindex',
            '-input {}'.format(db_fpath),
            '-iformat blastdb'
        ]
    )

    return makembindex_cmd
# end def



def create_reference_database(kromsatel_args):

    db_dirpath = os.path.join(kromsatel_args['outdir'], 'blast_database')
    fs.create_dir(db_dirpath)
    db_fpath = os.path.join(db_dirpath, 'kromsatel_blast_database')

    print('Creating the reference database for BLAST:\n  `{}`...'.format(db_fpath))
    _make_blast_db(kromsatel_args['reference_fpath'], db_fpath)
    print('Database: created')

    if kromsatel_args['use_index']:
        print('Indexing the database...')
        _index_database(db_fpath)
        print('Indexing: done')
    else:
        print('Index will not be created for the database.')
    # end if

    return db_fpath
# end def


def _make_blast_db(reference_fpath, db_fpath):
    makeblastdb_cmd = _configure_makeblastdb_cmd(reference_fpath, db_fpath)

    pipe = sp.Popen(makeblastdb_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print_err('\nCannot create blast database')
        print_err(stdout_stderr[1].decode('utf-8'))
        print_err(' Command: `{}`'.format(makeblastdb_cmd))
        platf_depend_exit(1)
    # end if
# end def


def _index_database(db_fpath):
    makembindex_cmd = _configure_makembindex_cmd(db_fpath)

    pipe = sp.Popen(makembindex_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print_err('\nCannot index the blast database `{}`'.format(db_fpath))
        print_err(stdout_stderr[1].decode('utf-8'))
        print_err(' Command: `{}`'.format(makembindex_cmd))
        platf_depend_exit(1)
    # end if
# end def


def _configure_blastn_cmd(query_fpath, db_fpath, blast_task, alignment_fpath):

    outfmt = 15 # Single-file BLAST JSON

    blast_cmd = ' '.join(
        [
            'blastn',
            '-query {}'.format(query_fpath),
            '-db {}'.format(db_fpath),
            '-task {}'.format(blast_task),
            '-evalue 1e-3',
            '-max_hsps 1', '-max_target_seqs 1',
            '-outfmt {}'.format(outfmt),
            '> {}'.format(alignment_fpath),
        ]
    )

    return blast_cmd
# end def


def blast_align(reads_chunk, kromsatel_args):

    query_fpath = os.path.join(
        kromsatel_args['outdir'],
        'kromsatel_query_{}.fasta'.format(os.getpid())
    )

    src.fastq.write_fastq2fasta(reads_chunk, query_fpath)

    alignment_fpath = os.path.join(
        kromsatel_args['outdir'],
        'kromsatel_alignment_{}.json'.format(os.getpid())
    )

    blast_cmd = _configure_blastn_cmd(
        query_fpath,
        kromsatel_args['db_fpath'],
        kromsatel_args['blast_task'],
        alignment_fpath
    )

    # Launch blastn
    pipe = sp.Popen(blast_cmd, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print_err('\nError while aligning a sequence against the database')
        print_err(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(1)
    # end if

    # fs.rm_temp_file(query_fpath)

    with open(alignment_fpath, 'rt') as alignment_file:
        aligmnents = json.load(alignment_file)
    # end with

    # fs.rm_temp_file(alignment_fpath)

    return aligmnents['BlastOutput2']
# end def blast_align
