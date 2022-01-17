
import os
import json
import subprocess as sp

import src.fastq
import src.filesystem as fs
from src.printing import getwt
from src.fatal_errors import FatalError


BLAST_TASKS = (
    'megablast',
    'dc-megablast',
    'blastn',
)

TASKS_SUPPORT_INDEXED_SEARCH = {
    BLAST_TASKS[0],
    BLAST_TASKS[2],
}


def get_blastplus_dependencies(krosatel_args):

    dependencies = [
        'blastn',
        'makeblastdb'
    ]

    if krosatel_args.use_index:
        dependencies.append('makembindex')
    # end if

    return dependencies
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
        error_msg = '\nError: program `{}` from BLAST+ toolkit is not installed.' \
            'If this error still occures although you have installed everything' \
            '  -- make sure that this program is added to PATH)' \
                .format(program)
        raise FatalError(error_msg)
    # end if
# end def check_blastn


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

    db_dirpath = os.path.join(kromsatel_args.outdir_path, 'blast_database')
    fs.create_dir(db_dirpath)
    db_fpath = os.path.join(db_dirpath, 'kromsatel_blast_database')

    print('{} - Creating the reference database for BLAST:\n  `{}`...'.format(getwt(), db_fpath))
    _make_blast_db(kromsatel_args.reference_fpath, db_fpath)
    print('{} - Database: created'.format(getwt()))

    if kromsatel_args.use_index:
        print('{} - Indexing the database...'.format(getwt()))
        _index_database(db_fpath)
        print('{} - Index: created'.format(getwt()))
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
        error_msg = '\nError: Cannot create blast database' \
            '{}\n Command: `{}`' \
                .format(stdout_stderr[1].decode('utf-8'), makeblastdb_cmd)
        raise FatalError(error_msg)
    # end if
# end def


def _index_database(db_fpath):
    makembindex_cmd = _configure_makembindex_cmd(db_fpath)

    pipe = sp.Popen(makembindex_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        error_msg = '\nError: Cannot index the blast database `{}`' \
            '{}\n Command: `{}`' \
                .format(db_fpath, stdout_stderr[1].decode('utf-8'), makembindex_cmd)
        raise FatalError(error_msg)
    # end if
# end def


def _configure_blastn_cmd_illumina(query_fpath, db_fpath, blast_task, use_index, alignment_fpath):

    outfmt = 15 # Single-file BLAST JSON

    if use_index:
        use_index_opt_value = 'true'
    else:
        use_index_opt_value = 'false'
    # end if

    blast_cmd = ' '.join(
        [
            'blastn',
            '-query {}'.format(query_fpath),
            '-db {}'.format(db_fpath),
            '-task {}'.format(blast_task),
            '-use_index {}'.format(use_index_opt_value),
            '-evalue 1e-3',
            '-gapopen 3 -gapextend 1',
            '-max_hsps 1', '-max_target_seqs 1',
            '-outfmt {}'.format(outfmt),
            '> {}'.format(alignment_fpath),
        ]
    )

    return blast_cmd
# end def


def _configure_blastn_cmd_nanopore(query_fpath, db_fpath, blast_task, use_index, alignment_fpath):

    outfmt = 15 # Single-file BLAST JSON

    if use_index:
        use_index_opt_value = 'true'
    else:
        use_index_opt_value = 'false'
    # end if

    blast_cmd = ' '.join(
        [
            'blastn',
            '-query {}'.format(query_fpath),
            '-db {}'.format(db_fpath),
            '-task {}'.format(blast_task),
            '-use_index {}'.format(use_index_opt_value),
            '-evalue 1e-3',
            '-max_hsps 3', '-max_target_seqs 1',
            '-outfmt {}'.format(outfmt),
            '> {}'.format(alignment_fpath),
        ]
    )

    return blast_cmd
# end def


def blast_align(reads_chunk, kromsatel_args):

    query_fpath = os.path.join(
        kromsatel_args.tmp_dir_path,
        'kromsatel_query_{}.fasta'.format(os.getpid())
    )

    src.fastq.write_fastq2fasta(reads_chunk, query_fpath)

    alignment_fpath = os.path.join(
        kromsatel_args.tmp_dir_path,
        'kromsatel_alignment_{}.json'.format(os.getpid())
    )

    if kromsatel_args.paired_mode:
        blast_cmd = _configure_blastn_cmd_illumina(
            query_fpath,
            kromsatel_args.db_fpath,
            kromsatel_args.blast_task,
            kromsatel_args.use_index,
            alignment_fpath
        )
    else:
        blast_cmd = _configure_blastn_cmd_nanopore(
            query_fpath,
            kromsatel_args.db_fpath,
            kromsatel_args.blast_task,
            kromsatel_args.use_index,
            alignment_fpath
        )
    # end if

    # Launch blastn
    pipe = sp.Popen(blast_cmd, shell=True, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        error_msg = '\nError: an error occured while performing BLAST search:' \
            '{}'.format(stdout_stderr[1].decode('utf-8'))
        raise FatalError(error_msg)
    # end if

    fs.rm_file_warn_on_error(query_fpath)

    with open(alignment_fpath, 'rt') as alignment_file:
        aligmnents = json.load(alignment_file)
    # end with

    fs.rm_file_warn_on_error(alignment_fpath)

    return aligmnents['BlastOutput2']
# end def blast_align
