## kromsatel

### Description

A script for preprocessing raw reads obtained using [ARTIC's protocol](https://artic.network/ncov-2019) for sequencing SARS-CoV-2 genome. Here, "preprocessing" stands for splitting concatemer reads into consistent fragments according to primer scheme described in the [protocol](https://artic.network/ncov-2019).

Brief description of the algorithm:

1. Align given read against amplicons using `blastn` program from BLAST+ toolkit (discontiguous megablast is used).

2. Extract alignments, which simultaneously:
    - Do start and/or end at 5'- and/or 3'-terminus of any amplicon.
    - Do not overlap within the read.

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

To create database of SARS-Cov-2 amplicons, do the following:

```
cd amplicons-db
makeblastdb -in nCoV-2019_amplicons.fasta -dbtype nucl -parse_seqids
```

When the database is created, you can proceed with cleaning.

#### Running kromsatel

```
./kromsatel.py <YOUR_READS> -d amplicons-db/nCoV-2019_amplicons.fasta 
```

You can also pass number of threads with `-t` option to align your reads in parallel. Further "cleaning" is still not parallelized :(.


#### Example
```
./kromsatel.py corona_reads.fastq -d amplicons-db/nCoV-2019_amplicons.fasta 
```

### Output file

For the example above, output file will have name `corona_reads_cleaned.fastq` and will be located is the same directory as `corona_reads.fastq`.

#### Read names

In output file, reads are named in following way:

```
  @<original_read_name>_<QSTART>_<QEND>
```

QSTART and QEND are 1-based coordinates of, correspondingly, start and end of an aligned fragment, which yields the output read (see [Description](#Description) section above for details).

Example of read name in an output file:
```
  @98786bd2-88a6-43ca-8c69-704992ad69cb_28-98
```

### Temporary files

While working, "kromsatel" creates temporary files and writes queries for "blastn" to them. If you terminate "kromsatel" or if it exits with error, they might be left after run.

Under Unix-like systems, these files are stored in `/tmp`, and under Windows -- in the working directory.

Naming schema of these files is following:

```
kromsatel_query_<PID>.fasta
```
, where PID is process ID of "kromsatel" process.
