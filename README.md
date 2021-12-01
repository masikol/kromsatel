## kromsatel

Current version is `1.6.b` (2021-12-01 edition).

### Description

A script for preprocessing raw reads obtained using [ARTIC's protocol](https://artic.network/ncov-2019) for sequencing SARS-CoV-2 genome. Here, "preprocessing" stands for splitting chimeric reads into consistent fragments according to primer scheme described in the [protocol](https://artic.network/ncov-2019) (or according to your own primer scheme).

Brief description of the algorithm:

1. Align given read against amplicons using `blastn` program from BLAST+ toolkit (discontiguous megablast is used).

2. Extract alignments, which do not overlap within the read and are long enough (see section "Options" about this "long enough").

3. Split the read into these aligned non-overlapping fragments (major amplicons are preferred).

### Dependencies

1. **Python 3** (https://www.python.org/). The script is tested on Python interpreter version 3.8.5.

  The script is written in Python 3X and won't work on Python 2X.

2. **BLAST+** toolkit.

   It can be downloaded here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

   **Installation:**

   - Linux: download tarball, unpack it and add `bin/` directory from unpacked tree to `PATH` variable.
     (this should also work for macOS)

   "kromsatel" has been tested on Linux with BLAST+ version 2.10.1+.

3. Following Python packages:

   - `numpy` (tested on version 1.19.2)
   - `pandas` (tested on version 1.1.3)

   They can be installed with following command: `pip3 install numpy pandas`

### Usage

#### Get started -- create the database

To create database of amplicons, run script `db-scripts/make-db.sh` in the following way:

```
bash db-scripts/make-db.sh <amplicons_fasta> <output_database_dir>

```
For example:
```
bash db-scripts/make-db.sh nCoV-2019_amplicons.fasta nCoV-2019_database
```

You can use other set of amplicons, for example the alternative one: `nCoV-2019-alt_amplicons.fasta` (it is generated using [this](https://github.com/ItokawaK/Alt_nCov2019_primers) primer set). You can also create your own amplicons for your primers using script `make-amplicons.sh` (see section "Creating your own amplicons" below).

When the database is created, you can proceed with cleaning.

#### Run kromsatel

#### Basic usage

```
./kromsatel.py <YOUR_READS> -d <DATABASE_PATH>
```

#### Options

```
-h (--help) -- print help message and exit.

-v (--version) -- print version and exit.

-o (--outdir) -- output directory.
    Default value is ./kromsatel_output

-d (--db) -- path to BLAST database of amplicons.
    Mandatory option (it sounds like an oxymoron but anyway).

-p (--primers-to-rm) -- CSV file containing primers and their names.
    It must be specified if you intend to remove primer sequences from reads.
    See section "Removing primer sequences" for details.

-t (--threads) -- number of threads to launch.
    Default: 1 thread.

--am -- minimum length of alignment against major amplicon.
    Shorter alignments will not be considered.
    Default: 100 bp.

--im -- minimum length of alignment against minor amplicon.
    Shorter alignments will not be considered.
    Default: 25 bp.

-c (--chunk-size) -- number of reads to blast within a single query.
    Default: 1000 reads.
```

#### Example
```
./kromsatel.py corona_reads.fastq \
  -d nCoV-2019_database/nCoV-2019_amplicons \
  -p primers/nCov-2019_primers.csv \
  -t 4 --am 150
```

And two commands together, for clarity:

Create database:

`bash db-scripts/make-db.sh nCoV-2019_amplicons.fasta nCoV-2019_database`

Run kromsatel

```
./kromsatel.py corona_reads.fastq \
  -d nCoV-2019_database/nCoV-2019_amplicons \
  -p primers/nCov-2019_primers.csv
```

You can pass multiple fastq files to kromsatel, like this:

```
./kromsatel.py \
  corona_reads_1.fastq corona_reads_2.fastq \
  -d nCoV-2019_database/nCoV-2019_amplicons \
  -p primers/nCov-2019_primers.csv
```

Or you can use wildcards (this will process all `.fastq.gz` file in directory `reads_dir`):

```
./kromsatel.py reads_dir/*.fastq.gz \
  -d nCoV-2019_database/nCoV-2019_amplicons \
  -p primers/nCov-2019_primers.csv
```

### Removing primer sequences

To make kromsatel remove sequences of primers, you must do the following:

1. When you create BLAST database, use fasta file (one among `amplicons-fasta/*.fasta`) WITHOUT words "with-primers" in it's name.

2. When you run kromsatel, pass CSV file containing primers with option `-p`.

And if you intend to keep primers sequence in reads -- do all the opposite: use "...with-primers" file of amplicons and omit `-p` option.

Not a convenient rule, but it is summarized in the table below:

|                                | `-p` specified | no `-p` option |
| ------------------------------ |:--------------:| --------------:|
| `amplcions.fasta`              |       OK       |   NOT CORRECT  |
| `amplcions_with-primers.fasta` |   NOT CORRECT  |       OK       |


### Output file

For the example above, output file will have name `corona_reads_cleaned.fastq` and will be located is the same directory as `corona_reads.fastq`.

#### Read names

In output file, reads are named in following way:

```
  @<original_read_name>_<QSTART>-<QEND>
```

QSTART and QEND are 1-based coordinates of, correspondingly, start and end of an aligned fragment, which yields the output read (see [Description](#Description) section above for details).

Example of read name in an output file:
```
  @98786bd2-88a6-43ca-8c69-704992ad69cb_28-98
```

### Creating your own amplicons

You can create your own amplicons for your primers using script `db-scripts/make-amplicons.sh`.

To do this, you should configure proper CSV (`primer_name,primer_seq`) file containing your primers. You can use file `primers/nCov-2019_primers.csv` as example. "Proper" means the following:

1. Primers names and sequences must be separated with comma.

2. Order of primers is crucial, so keep order exactly as in this example file.

#### Dependencies

This script depends on [seqkit](https://github.com/shenwei356/seqkit).

#### Input:

1. CSV File containing primers (see file `primers/nCov-2019_primers.csv` for example);

2. Genome (in fasta format) of damned coronavirus (I use this one: [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512)).

#### Output:

Script "make-amplicons.sh" outputs two fasta files containing amplicons:

1. File containing amplicons whith primers sequences removed, e.g. `nCoV-2019_amplicons.fasta`.

2. File containing amplicons where primers sequences are not removed, e.g. `nCoV-2019_amplicons_with-primers.fasta`.

Any of these output files are ready for database creation, just like `amplicons-fasta/nCoV-2019_amplicons.fasta`.

#### Options:
```
 -h -- print help message and exit;
 -p -- CSV file with primers  (mandatory);
 -g -- fasta file with genome (mandatory);
 -o -- prefix of output fasta files;
       default value: ./my-amplicons
 -s -- path to seqkit executable
       (if seqkit is in your PATH, just omit this option);
 -i -- minimum length of a minor amplicon to output;
       default value: 25.
 ```

#### Usage:
```
  bash make-amplicons.sh <-p primers_csv> <-g genome_fasta> [-o STR] [-s STR] [-i INT]
```
#### Examples:

Case 1. `seqkit` is in PATH
```
  bash make-amplicons.sh -p my_primers.csv -g Wuhan-Hu-1-compele-genome.fasta \
  -o amplicons-fasta/my_amplicons.fasta
```
Case 2. `seqkit` is not in PATH
```
  bash make-amplicons.sh -p my_primers.csv -g Wuhan-Hu-1-compele-genome.fasta \
  -o amplicons-fasta/my_amplicons.fasta -s /home/me/seqkit/bin/seqkit
```

### Kromsatel temporary files

While working, "kromsatel" creates temporary files and writes queries for "blastn" to them. If you terminate "kromsatel" or if it exits with error, they might be left after run.

Under Unix-like systems, these files are stored in `/tmp`, and under Windows -- in user's temporary directory.

Naming schema of these files is following:

```
kromsatel_query_<PID>.fasta
```
, where PID is process ID of "kromsatel" process.

