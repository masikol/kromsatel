# kromsatel

Current version is `1.7.a_dev` (2021-XX-XX edition).

## Description (TODO: depracated)

Deprecation notice: now kromsatel can process only Illumina paired-end reads.

Kromsatel is a script for preprocessing raw reads obtained using [ARTIC's protocol](https://artic.network/ncov-2019) for sequencing SARS-CoV-2 genome. Here, "preprocessing" stands for splitting chimeric reads into consistent fragments according to primer scheme described in the [protocol](https://artic.network/ncov-2019) (or according to your own primer scheme).

Brief description of the algorithm:

1. Align given read against amplicons using `blastn` program from BLAST+ toolkit (discontiguous megablast is used).

2. Extract alignments, which do not overlap within the read and are long enough (see section "Options" about this "long enough").

3. Split the read into these aligned non-overlapping fragments (major amplicons are preferred).

### Dependencies

1. **Python 3** (https://www.python.org/). The script is tested on Python interpreter version 3.8.5.

  The script is written in Python 3 and won't work on Python 2.

2. **BLAST+** toolkit.

   It can be downloaded [here](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

   BLAST+ Installation:

   - Linux: download tarball, unpack it and add `bin/` directory from unpacked tree to `PATH` variable.
     (this should also work for macOS)

   Kromsatel has been tested on Linux with BLAST+ version 2.12.0+.

## Usage

### Arguments

Mandatory arguments are marked with `*`.

```
  -h (--help) -- print help message and exit.

  -v (--version) -- print version and exit.

Input:

* -1 (--reads-R1) -- a fastq file of forward reads.
      The file may be gzipped.

* -2 (--reads-R2) -- a fastq file of reverse reads.
      The file may be gzipped.

  -u (--reads-unpaired) -- a fastq file of unpaired reads.
      The file may be gzipped.
      (TODO: Not supported yet)

* -p (--primers) -- a CSV file of primer names and sequences.
      This must be two-columns CSV file without header.

* -r (--reference) -- a fasta file of reference sequence.

Output:

-o (--outdir) -- output directory.
    Default value is './kromsatel_output'

Miscellaneous:

-t (--threads) -- number of threads to launch.
    Default: 1 thread.

-m (--min-len) -- minimum length of an output read.
    Default: 25 bp.

-a (--blast-task) -- BLASTn task to launch.
    Allowed values: 'megablast', 'dc-megablast', 'blastn'.
    Default is 'megablast'.

-c (--chunk-size) -- number of reads to blast within a single query.
    The larger is the chunk size, the higher is the memory consumption.
    Default: 1000 reads.
```

### Examples
```
./kromsatel.py \
    -1 20_S30_L001_R1_001.fastq.gz \
    -2 20_S30_L001_R2_001.fastq.gz \
    -p primers/nCov-2019_primers.csv \
    -r reference/Wuhan-Hu-1-compele-genome.fasta \
    -o 20_s30_outdir
```

## Output files

TODO

### Read names (TODO: deprecated)

In output file, reads are named in following way:

```
  @<original_read_name>_<QSTART>-<QEND>
```

QSTART and QEND are 1-based coordinates of, correspondingly, start and end of an aligned fragment, which yields the output read (see [Description](#Description) section above for details).

Example of read name in an output file:
```
  @98786bd2-88a6-43ca-8c69-704992ad69cb_28-98
```
