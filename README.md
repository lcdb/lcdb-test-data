# LCDB test data

This repository stores small amounts of example data that can be used for
testing along with the code to reproduce the data.

## Building data sets

See `lcdb_test_data/Snakefile` for details. To regenerate data, install this
package and run `build_example_data`. This will build a conda environment with
the required dependencies and use that environment to run the Snakefile.


```bash
python setup.py develop
build_example_data -h
build_example_data data /path/to/example/data -j8
```

## RNA-seq data

RNA-seq test data currently consists of four, 48-bp PE samples from GEO
accession GSE49587.  This was an RNA-seq experiment WT and Smn mutant
Drosophila larvae.

In addition to the small FASTQs, we also have matching genome, transcriptome,
and GTF files.

label   | SRA accession | treatment         |
--------|---------------|-------------------|
sample1 | SRR948304     | WT rep 1          |
sample2 | SRR948305     | WT rep 2          |
sample3 | SRR948306     | Smn mutant, rep 1 |
sample4 | SRR948307     | Smn mutant, rep 2 |


## ChIP-seq data

ChIP-seq data consists of three IP/input pairs from GSE38594.

label  | SRA accession | description                 |
-------|---------------------------------------------|
input1 | SRR504958     | wing disc input, rep1       |
input2 | SRR504959     | wing disc input, rep2       |
input3 | SRR504959     | embyro input, rep1          |
ip1    | SRR504955     | GAF ChIP in wing disc, rep1 |
ip2    | SRR504956     | GAF ChIP in wing disc, rep2 |
ip3    | SRR504946     | GAF ChIP in embyro, rep1    |
ip4    | SRR504947     | GAF ChIP in embyro, rep2    |


