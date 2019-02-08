# egfr-v3-determiner #

## What is does ##

Estimates the number of EGFR-vIII and EGFR-wt reads in a BAM file:

```
$ bin/egfr-v3-determiner -r hg38 tmp/test_001.bam
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
  -r, --reference-build [hg19|hg38]
                                  Used reference genome (needed for EGFR exon
                                  coordinates)  [required]
  --help                          Show this message and exit.
```
