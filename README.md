# kromsatel

Current version is `1.7.d_dev` (2022-01-XX edition).

Kromsatel is a program which preprocesses ("cleans") raw reads of amplicon sequencing.

For example, when the genome of SARS-CoV-2 is sequenced, reads of amplicons produced using [ARTIC](https://artic.network/ncov-2019) protocol can be processed with kromsatel prior to downstream genome analysis. Any other amplicon protocol is (most likely) also acceptable; however, only ARTIC was tested.

## Description

Kromsatel can process Illumina (paired-end or single-end) and Nanopore read sets.

### Preprocessing: specifics

For kromsatel, the "preprocessing" stands for the following:

1) Remove primer sequences from reads.

2) Trim reads so that any unaligned part of read is removed.

3) Discriminate reads coming from major, minor, and non-specific amplicons (see the picture below, TODO).

4) Split chimeric reads into consistent fragments (only Nanopore reads).

## Dependencies

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
  -h (--help) -- print the help message and exit.

  -v (--version) -- print version and exit.

Input:

  -1 (--reads-R1) -- a fastq file of forward reads.
      The file may be gzipped.

  -2 (--reads-R2) -- a fastq file of reverse reads.
      The file may be gzipped.

  -u (--reads-unpaired) -- a fastq file of unpaired reads.
      The file may be gzipped.

* -p (--primers) -- a CSV file of primer names and sequences.
      This file must be a two-column CSV file without header.

* -r (--reference) -- a fasta file of reference sequence.

Output:

  -o (--outdir) -- output directory.
      Default value is './kromsatel_output'

  -m (--min-len) -- minimum length of an output read.
      Default: 25 bp.

Computational resources:

  -t (--threads) -- number of threads to launch.
      Default: 1 thread.

Advanced:

  -k (--blast-task) -- BLASTn task to launch.
      Allowed values: 'megablast', 'dc-megablast', 'blastn'.
      Default is 'megablast'.

  -c (--chunk-size) -- number of reads to blast within a single query.
      The larger is the chunk size, the higher is the memory consumption.
      Default: 1000 reads.

  --crop-len -- number of nucleotides to crop from end of reads
      originating from a non-specific amplicon.
      Default: 'auto' (maximum primer length).

  --primer-5ext -- size of 5' primer coordinates extention.
      Low (< 2 bp) values of this parameter may result in extra major alignments
      classified as uncertain, and vice versa for high values (> 10 bp).
      Default: 5 bp.

  --use-index -- Whether to use BLAST index.
      Permitted values: auto, true, false.
      "auto" mode: true for megablast and blastn, false for dc-megablast.
      Default: false.
```

### Examples

#### Short (e.g. Illumina) paired-end reads

```
./kromsatel.py \
    -1 20_S30_L001_R1_001.fastq.gz \
    -2 20_S30_L001_R2_001.fastq.gz \
    -p primers/nCov-2019_primers.csv \
    -r reference/Wuhan-Hu-1-compele-genome.fasta \
    -o 20_S30_outdir
```

#### Short (e.g. Illumina) single-end reads
```
./kromsatel.py \
    -1 20_S30_L001_R1_001.fastq.gz \
    -p primers/nCov-2019_primers.csv \
    -r reference/Wuhan-Hu-1-compele-genome.fasta \
    -o 20_S30_outdir
```

#### Long (e.g. nanopore) reads
```
./kromsatel.py \
    -u Wuhan-Hu-1.fastq.gz \
    -p primers/nCov-2019_primers.csv \
    -r reference/Wuhan-Hu-1-compele-genome.fasta \
    -o Wuhan-Hu-1_outdir
```

#### With additional options
```
./kromsatel.py \
    -1 20_S30_L001_R1_001.fastq.gz \
    -2 20_S30_L001_R2_001.fastq.gz \
    -p primers/nCov-2019_primers.csv \
    -r reference/Wuhan-Hu-1-compele-genome.fasta \
    -o 20_S30_outdir \
    -k dc-megablast \
    -c 2000 \
    -m 50 \
    -t 4 \
    --crop-len 27 \
    --primer-5ext 3
```

## Output files

#### Paired-end reads

Let input read files be `20_S30_L001_R1_001.fastq.gz` and `20_S30_L001_R2_001.fastq.gz`. For them, kromsatel will produce the following output files:

Pair of files of reads coming from major amplicons:
```
20_S30_L001_R1_001_major.fastq.gz
20_S30_L001_R2_001_major.fastq.gz
```

Pair of files of reads coming from minor amplicons:
```
20_S30_L001_R1_001_minor.fastq.gz
20_S30_L001_R2_001_minor.fastq.gz
```

Pair of files of reads coming from non-specific or undetected ampicons:
```
20_S30_L001_R1_001_uncertain.fastq.gz
20_S30_L001_R2_001_uncertain.fastq.gz
```

Pair of files of reads which lost their pair during preprocessing:
```
20_S30_L001_R1_001_unpaired.fastq.gz
20_S30_L001_R2_001_unpaired.fastq.gz
```
#### Single-end reads

Let input read file be `all_pass_15.fastq.gz`. For them, kromsatel will produce the following output files:

File if reads coming from major amplicons:

```
all_pass_15_major.fastq.gz
```

File if reads coming from minor amplicons:

```
all_pass_15_minor.fastq.gz
```

File if of reads coming from non-specific or undetected ampicons:

```
all_pass_15_uncertain.fastq.gz
```

### Read names

#### Illumina (paired-end) reads

Kromsatel keeps read headers unchanged.

#### Nanopore reads

For Nanopore data, kromsatel modifies headers of output reads. Thus, the header of an output read will looklike this:

```
  @<ORIGINAL_READ_NAME>_<QSTART>-<QEND> [runid,sampleid,etc...]
```

, where QSTART and QEND are 1-based coordinates of, correspondingly, start and end of the aligned fragment, which yields the output read.

An example of a modified read name in an output file:
```
  @98786bd2-88a6-43ca-8c69-704992ad69cb_28-98 runid=...,sampleid=...,etc...
```
