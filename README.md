# egfr-v3-determiner #

## What is does ##

Estimates the number of EGFR-vIII and EGFR-wt reads in a BAM file:

```
$ egfr-v3-determiner -r hg38 tmp/test_001.bam
```

Will result in a text file like this:

| sample | wt-reads | vIII-reads |
|--------|----------|------------|
| tmp/test_001.bam | 0 | 1 |

## Installation ##

```
git clone https://github.com/yhoogstrate/egfr-v3-determiner.git
cd egfr-v3-determiner

virtualenv -p python3 .venv
source .venv/bin/activate

python setup.py install

egfr-v3-determiner --help
```

## Usage ##

```
Usage: egfr-v3-determiner [OPTIONS] [INPUT_BAM]...

Options:
  --version                       Show the version and exit.
  -r, --reference-build [hg19|hg38]
                                  Used reference genome (needed for EGFR exon
                                  coordinates)  [required]
  -s, --spliced-reads-only        If paired end reads with an insert size
                                  longer than 801 bases can be expected, wild-
                                  type exon-1 to exon-8 covering reads can
                                  be expected. Enabling this flag only uses
                                  spliced reads for vIII determination.
  -n, --read-names                Report all read-names instead of the read
                                  counts.
  -i, --include-interchromosomal  Include paired-end reads that have an
                                  interchromosomal mapped mate (disabled by
                                  default).
  -d, --dataset-suffix TEXT       Adds this suffix to the column names; tabs
                                  and newlines not allowed.
  --help                          Show this message and exit.
```
