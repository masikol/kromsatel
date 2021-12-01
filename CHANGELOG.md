# kromsatel changelog

## 2021-12-01 edition

Add option `-o/--outdir` (output directory) to kromsatel.py.

Version change: `1.5.a --> 1.6.a`.

## 2021-05-21 edition

Improved algorithm for removing primer sequences. Although it slows kromsatel down twofold, it removes primers pretty well.

See section ["Removing primer sequences"](https://github.com/masikol/kromsatel/blob/main/README.md#removing-primer-sequences) in README for details.

Version change: `1.4.a --> 1.5.a`.

## 2021-03-29 edition

Added possibility for kromsatel to remove primers sequences from reads being processed. To do that, you should create your BLAST database using fasta file of amplicons, in which primer sequences are removed. Those files do not contain phrase `_with-primers` in their names (see folder `amplicons-fasta/` in the repo).

Added script `db-scripts/make-db.sh` version `1.0.a`

Version change: `db-scripts/make-amplicons.sh`: `1.0.a --> 1.1.a`.

## 2021-03-05 edition

- Improved performance ~ 2.3 times.
- Enabled parallel read processing (`-t`).
- Changed the very algorithm: it no longer looks if an alignment spans from primer to primer (or "starts" at one primer and interrupts somewhere in beteen). Now it just removes too short alignments (see options `--am` and `--im`).
- Added options `-c`, `--am`, `--im`.

Version change: `1.3.b --> 1.4.a`.

## 2021-01-27 afternoon edition

Made the algorithm less brutal to those reads which do not span from primer to primer. Now kromsatel allows offset of 10 bp from both ends (this offset is no longer actual).

Version change: `1.3.a --> 1.3.b`.

## 2021-01-27 edition

Added script `make-amplicons.sh` and new (alternative) set of amplicons: `amplicons-db/nCoV-2019-alt_amplicons.fasta`.

Primers for this new amplicons are from [here](https://github.com/ItokawaK/Alt_nCov2019_primers).

## 2021-01-24

Moved kromsatel from [cager-misc](https://github.com/masikol/cager-misc) to this separate repo.

## 2021-01-06 edition

- kromsatel: fixed bug that would cause kromsatel to prefer minor amplicons over major ones.

### Version changes:

- kromsatel: `1.2.e -> 1.3.a`

## 2020-12-02 edition

- kromsatel: output files naming changed: now "cleaned" is suffix, not prefix.

### Version changes:

- kromsatel: `1.2.d -> 1.2.e`

## 2020-11-14 edition

- kromsatel now works 2 times faster.

### Version changes:

- kromsatel: `1.2.c -> 1.2.d`

## 2020-11-13 edition

- kromsatel: performance improved.

### Version changes:

- kromsatel: `1.2.b -> 1.2.c`

## 2020-11-12 edition

- kromsatel version `1.2.a` added, then added version `1.2.b` with code commented and some optimized operations.
