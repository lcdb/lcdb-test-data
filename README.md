# LCDB test data

This repository stores small amounts of example data that can be used for
testing along with the code to reproduce the data.

It currently consists of four, 48-bp PE samples from GEO accession GSE49587.
This was an RNA-seq experiment WT and Smn mutant Drosophila larvae. These
FASTQs were first mapped to the full dm6 assembly. Then a subset of reads is taken
from the first megabase of chromosomes 2L and 2R. These are added to the final
subset FASTQ. An additional 10000 or so unmapped reads are also added.

In addition to the small FASTQs, we also have matching genome, transcriptome,
and GTF files.


label   | SRA accession | treatment         |
--------|---------------|-------------------|
sample1 | SRR948304     | WT rep 1          |
sample2 | SRR948305     | WT rep 2          |
sample3 | SRR948306     | Smn mutant, rep 1 |
sample4 | SRR948307     | Smn mutant, rep 2 |


See `lcdb_test_data/Snakefile` for details. To regenerate data, install this
package and run `build_example_data`. This will build a conda environment with
the required dependencies and use that environment to run the Snakefile.


```
python setup.py develop
build_example_data -h
build_example_data data /path/to/example/data -j8
```
