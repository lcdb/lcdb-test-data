# LCDB test data

This repository stores small amounts of example data that can be used for
testing along with the code to reproduce the data.

## Building data sets

The Snakefile assumes an environment called `lcdb-test-data`, created like this
(assuming bioconda channel is set up):

```bash
conda create -n lcdb-test-data --file requirements.txt
```

Run the Snakefile as normal. Note that this will be downloading and creating
10s of GB of data, so you will want to specify the working directory as
somewhere with a fair amount of room (use the `--directory` or `-d` arg for
Snakemake).

You can use the `WRAPPER_SLURM` file for submitting to a SLURM cluster. E.g.,
on NIH's biowulf:

```
sbatch WRAPPER_SLURM Snakefile -d /path/to/full/data
```

When complete, run the `cp-data-to-repo.sh` script to just copy over the small
example data to this repo, and then make the commits.

## Strategy

The full reference genome, annotations, and transcriptome are downloaded.
FASTQs are aligned to the full reference. The resulting BAMs are parsed and
downsampled across the region indicated in the `LIMITS.bed` file (created by the
workflow; configured in the `limits` rule) in such a way that read pairs are
handled correctly.

The reference genome, transcriptome, and annotations are similarly subset to
the region indicated in `LIMITS.bed`.

To avoid too-sparse data, the BAMs are first subset by the restricted region
and then downsampled. The amount of downsampling and ratio of
mapped-to-unmapped reads are set by the `mapped_n_config` and
`unmapped_n_config` dictionaries There are "small" and "tiny" versions of each.

These newly-subset FASTQs are then re-aligned to the genome

In some cases, you may find that tests do not not have enough reads. In this
case, reset the values in those dictionaries and then force-rerun the
`chipseq_small_fastq` and `rnaseq_small_fastq` rules:

```
snakemake -d /path/to/output --forcerun chipseq_small_fastq rnaseq_small_fastq
```

A slightly more extreme adjustment would be use change the coordinates in the
`limits` rule, and then force-run that rule.

```
snakemake -d /path/to/output --forcerun limits
```

## RNA-seq data

RNA-seq test data currently consists of four, 48-bp PE samples from GEO
accession GSE49587.  This was an RNA-seq experiment WT and Smn mutant
Drosophila larvae.


label   | SRA accession | treatment         |
--------|---------------|-------------------|
sample1 | SRR948304     | WT rep 1          |
sample2 | SRR948305     | WT rep 2          |
sample3 | SRR948306     | Smn mutant, rep 1 |
sample4 | SRR948307     | Smn mutant, rep 2 |


## ChIP-seq data

ChIP-seq data consists of three IP/input pairs from GSE38594.

label  | SRA accession | description                 |
-------|---------------|-----------------------------|
input1 | SRR504958     | wing disc input, rep1       |
input2 | SRR504959     | wing disc input, rep2       |
input3 | SRR504959     | embyro input, rep1          |
ip1    | SRR504955     | GAF ChIP in wing disc, rep1 |
ip2    | SRR504956     | GAF ChIP in wing disc, rep2 |
ip3    | SRR504946     | GAF ChIP in embyro, rep1    |
ip4    | SRR504947     | GAF ChIP in embyro, rep2    |

