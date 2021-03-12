# egfr-v3-determiner #

EGFRvIII is the most common mutation found in the EGFR gene in glioblastomma.
It originates from a genomic deletion or multiple inversion resulting in the
loss of exons 2 - 7 at RNA level.

Using RNA-seq, EGFRvIII can be determined based on reads that:

 - (1) fall exactly over the splice junction (use `-s` / `--spliced-only`)
 - (2) are 'spanning' the splice junction and are thus mapped perfectly within exons 1 and 8

**The only thing that you need are properly aligned RNA-seq BAM files and this tool.**

Please notice that type 2 reads can also arise with other structural variants
in combination with large insert sizes. For instance, if exons 2 - 6 are
deleted and the length of the inner sequence of the RNA frament is longer than
exon 7, it could be wrongly interpreted as EGFRvIII read.


## What is does ##

We designed a small python tool for estimnating the read counts and/or extracting
the actual read names that allows to further analysed the sequencing data.

Estimates the number of EGFR-vIII and EGFR-wt reads from a BAM file directly:

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

### Genomic reference ###

The genomic reference (hg19/hg39) can be changed using the `-r` or
`--reference-build` argument. Genomic references starting with '>1' rather
than 'chr1' will be automatically resolved.

### Interchromosomal reads ###

It may happen that one of the mates splices perfectly over exons 1 - 8 but
that it's mate is mapped to another chromosome. As these are likely derived
from other stuctural variants or odd reads, they are by default excluded.

It may nevertheless, be interesting to analyse these reads, for instance if
there is a suspicion for other structural variants. By using the `-i` /
`--include-interchromosomal` argument, these reads will be included.

### Suffix output column names ###

It may be convenient to add a suffix to the column name in the output, for
instance '-EGFRvIII-reads'. This can be achieved by using the `-d` / 
`--dataset-suffix` argument:

```
egfr-v3-determiner -d '-EGFRvIII-reads' rna-seq-sample-1.bam rna-seq-sample-2.bam > EGFRvIII.counts.txt
```

## I NEED HELP / I FOUND A BUG / I WANT A FEATURE ##

Great, please do so :) I am more than happy to help.
If you find some odd reads that I need to take a look at, feel free
to send them to me.

If you prefer private communication you can e-mail me at:

y {dot} hoogstrate {at} erasmusmc {.} nl


## LICENSE ##

This is FREE software without liability or warranty, following the GPL-3
software license.


