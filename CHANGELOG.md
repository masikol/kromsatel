# kromsatel changelog

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