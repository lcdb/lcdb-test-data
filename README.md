# LCDB test data

This repository stores small amounts of example data that can be used for
testing along with the code to reproduce the data.

It currently consists of four, 48-bp PE samples from GEO accession GSE49587.
This was an RNA-seq experiment WT and Smn mutant Drosophila larvae. These
FASTQs were mapped to the first megabase of chromosomes 2L and 2R from the dm6
assembly. All mapped reads and 10% of unmapped reads are then extracted from
the BAMs to make new, minimal FASTQ files.

See `lcdb_test_data/Snakefile` for details. To regenerate data, install this
package and run `build_example_data`. This will build a conda environment with
the required dependencies and use that environment to run the Snakefile.

```
python setup.py develop
build_example_data -h
build_example_data data /path/to/example/data -j8
```
