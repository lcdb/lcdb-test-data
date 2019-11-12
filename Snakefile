import gzip
import os
from Bio import SeqIO
from Bio.Seq import Seq
from textwrap import dedent

# required to avoid near-simultaneous timestamps that confuse snakemake
shell.prefix('sleep 2; source activate env/; ')

# This is required for running sratoolkit on biowulf; if you don't need it then
# set to empty string
VDB_CONFIG_PRELUDE = 'export VDB_CONFIG=/usr/local/apps/ncbi/config/biowulf.kfg'

rnaseq_accessions = {
    'sample1': 'SRR948304',
    'sample2': 'SRR948305',
    'sample3': 'SRR948306',
    'sample4': 'SRR948307',
}

chipseq_accessions = {
    'input_1': 'SRR504958',
    'input_2': 'SRR504959',
    'input_3': 'SRR504954',
    'ip_1': 'SRR504955',
    'ip_2': 'SRR504956',
    'ip_3': 'SRR504946',
    'ip_4': 'SRR504948',
}

cutnrun_accessions = {
    'cutnrun_ip_1': 'SRR7988793',
    'cutnrun_ip_2': 'SRR7988784',
}

mapped_n_config = dict(small=2000000, tiny=10000)
unmapped_n_config = dict(small=1000, tiny=100)
multimapped_n_config = dict(small=5000, tiny=500)

n = range(1, 5)
rule all:
    input:
        [
            'data/annotation/dm6.small.refflat',
            'data/seq/dm6.small.fa',
            'data/seq/dm6.small.transcriptome.fa',
            'data/LIMIT.bed',
        ]
        + expand('data/rnaseq_samples/{sample}/{sample}.{size}_R{N}.fastq.gz',
                n=n, N=[1,2], size=['full', 'small', 'tiny'], sample=rnaseq_accessions.keys())
        + expand('data/rnaseq_samples/{sample}/{sample}.{size}.{r}.sorted.bam',
                size=['full', 'small', 'tiny'], n=n, r=['single', 'paired'], sample=rnaseq_accessions.keys())
        + expand('data/chipseq_samples/{sample}/{sample}.{size}_R1.fastq.gz',
                size=['full', 'small', 'tiny'], sample=chipseq_accessions.keys())
        + expand('data/chipseq_samples/{sample}/{sample}.{size}.single.sorted.bam',
                size=['full', 'small', 'tiny'], sample=chipseq_accessions.keys())
        + expand('data/cutnrun_samples/{sample}/{sample}.{size}_R{N}.fastq.gz',
                N=[1,2], size=['full', 'small', 'tiny'], sample=cutnrun_accessions.keys())
        + expand('data/cutnrun_samples/{sample}/{sample}.{size}.paired.sorted.bam',
                size=['full', 'small', 'tiny'], sample=cutnrun_accessions.keys())


# ----------------------------------------------------------------------------
# Create a BED file that will be used to subset GTF and FASTA files
rule limits:
    output: 'data/LIMIT.bed'
    shell:
        'echo "chr2L	0	1000000	chr2L" > {output}; '
        'echo "chr2R	0	1000000	chr2R" >> {output}'


# ----------------------------------------------------------------------------
# Download FlyBase GTF
rule prep_gtf:
    output: 'data/annotation/dm6.full.gtf'
    shell:
        'wget --no-clobber -q '
        '-O- '
        'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.11_FB2016_03/gtf/dmel-all-r6.11.gtf.gz > tmp.gtf.gz '
        '&& zcat tmp.gtf.gz '
        '| bedtools sort -i stdin '
        '| grep exon '
        """| awk '{{print "chr"$0}}'  > {output} """
        '&& rm tmp.gtf.gz '


# ----------------------------------------------------------------------------
# Subset GTF based on limits
rule prep_small_gtf:
    input:
        gtf=rules.prep_gtf.output,
        limit=rules.limits.output
    output: 'data/annotation/dm6.small.gtf'
    shell:
        'bedtools intersect -a {input.gtf} -b {input.limit} > {output} '


# ----------------------------------------------------------------------------
# Download Flybase transcriptome FASTA
rule prep_transcriptome:
    input: rules.prep_gtf.output
    output: 'data/seq/dm6.full.transcriptome.fa'
    shell:
        'wget --no-clobber -q '
        '-O- '
        'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.11_FB2016_03/fasta/dmel-all-transcript-r6.11.fasta.gz '
        '| gzip -d -c > {output} '

# ----------------------------------------------------------------------------
# Subset transcriptome based on transcript IDs retained in the subsetted GTF
rule prep_small_transcriptome:
    input:
        gtf=rules.prep_small_gtf.output,
        fasta=rules.prep_transcriptome.output
    output: 'data/seq/dm6.small.transcriptome.fa'
    run:
        from Bio import SeqIO
        import gffutils
        features = gffutils.iterators.DataIterator(str(input.gtf))
        keep = set([i.attributes['transcript_id'][0] for i in features])
        parser = SeqIO.parse(str(input.fasta), 'fasta')
        recs = []
        for rec in parser:
            if rec.name in keep:
                recs.append(rec)
        with open(output[0], 'w') as fout:
            SeqIO.write(recs, fout, 'fasta')

# ----------------------------------------------------------------------------
# Convert small GTF to refflat
rule gtftorefflat:
    input: rules.prep_small_gtf.output
    output: 'data/annotation/dm6.small.refflat'
    shell:
        'gtfToGenePred {input} {output}.tmp '
        '&& paste <(cut -f1 {output}.tmp) {output}.tmp > {output} '
        '&& rm {output}.tmp'


# ----------------------------------------------------------------------------
# Download full fasta
rule prep_fasta:
    input: rules.limits.output
    output: 'data/seq/dm6.full.fa'
    shell:
        'wget --no-clobber -q '
        '-O- ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.11_FB2016_03/fasta/dmel-all-chromosome-r6.11.fasta.gz '
        '| gunzip -c '
        '| sed "s/>/>chr/g" > {output} '

# ----------------------------------------------------------------------------
# Subset genome fasta
rule prep_small_fasta:
    input:
        fasta=rules.prep_fasta.output,
        limits=rules.limits.output
    output: 'data/seq/dm6.small.fa'
    shell:
        'bedtools getfasta -fi {input.fasta} -bed {input.limits} | '
        '''awk -F ":" '/^>/{{print $1; next}}{{print}}' > {output} '''


rule download_rnaseq_fastqs:
    output:
        fastq_R1='data/rnaseq_samples/{sample}/{sample}.full_R1.fastq.gz',
        fastq_R2='data/rnaseq_samples/{sample}/{sample}.full_R2.fastq.gz'
    run:
        accession = rnaseq_accessions[wildcards.sample]
        outdir = 'data/raw'
        shell('{VDB_CONFIG_PRELUDE}; fastq-dump --outdir {outdir} {accession} --split-files')
        shell('gzip -c {outdir}/{accession}_1.fastq > {output.fastq_R1}')
        shell('gzip -c {outdir}/{accession}_2.fastq > {output.fastq_R2}')

        # clean up raw fastqs
        shell('rm {outdir}/{accession}_1.fastq')
        shell('rm {outdir}/{accession}_2.fastq')

rule download_chipseq_fastqs:
    output:
        fastq_R1='data/chipseq_samples/{sample}/{sample}.full_R1.fastq.gz',
    run:
        accession = chipseq_accessions[wildcards.sample]
        outdir = 'data/raw'
        shell('{VDB_CONFIG_PRELUDE}; fastq-dump --outdir {outdir} {accession}')
        shell('gzip -c data/raw/{accession}.fastq > {output.fastq_R1}')

        # clean up raw fastqs
        shell('rm {outdir}/{accession}.fastq')


rule download_cutnrun_fastqs:
    output:
        fastq_R1='data/cutnrun_samples/{sample}/{sample}.full_R1.fastq.gz',
        fastq_R2='data/cutnrun_samples/{sample}/{sample}.full_R2.fastq.gz',
    run:
        accession = cutnrun_accessions[wildcards.sample]
        outdir = 'data/raw'
        shell('{VDB_CONFIG_PRELUDE}; fastq-dump --outdir {outdir} {accession} --split-files')
        shell('gzip -c {outdir}/{accession}_1.fastq > {output.fastq_R1}')
        shell('gzip -c {outdir}/{accession}_2.fastq > {output.fastq_R2}')

        # clean up raw fastqs
        shell('rm {outdir}/{accession}_1.fastq')
        shell('rm {outdir}/{accession}_2.fastq')


rule download_all_fastqs:
    input:
        expand('data/chipseq_samples/{sample}/{sample}.full_R1.fastq.gz', sample=chipseq_accessions.keys()),
        expand('data/rnaseq_samples/{sample}/{sample}.full_R{n}.fastq.gz', sample=rnaseq_accessions.keys(), n=[1,2]),
        expand('data/cutnrun_samples/{sample}/{sample}.full_R{n}.fastq.gz', sample=cutnrun_accessions.keys(), n=[1,2]),

# ----------------------------------------------------------------------------
# HISAT2 index
rule hisat_index:
    input: rules.prep_fasta.output
    output: expand('data/seq/dm6.full.{n}.ht2', n=range(1,8))
    params: index='data/seq/dm6.full'
    log: 'data/seq/dm6.small.ht2.log'
    shell:
        'hisat2-build {input} {params.index} &> {log}'

# ----------------------------------------------------------------------------
# Bowtie2 index
rule bowtie2_index:
    input: rules.prep_fasta.output
    output: expand('data/seq/dm6.full.{n}.bt2', n=range(1,2))
    params: index='data/seq/dm6.full'
    log: 'data/seq/dm6.small.bt2.log'
    shell:
        'bowtie2-build {input} {params.index} &> {log}'

# ----------------------------------------------------------------------------
# HISAT2 align.
#
# Note we're creating both SE and PE bams in serial rather than parallel
# (simplifies the snakefile)
rule hisat_align:
    input:
        index=expand('data/seq/dm6.full.{n}.ht2', n=range(1,8)),
        fastq_R1='data/rnaseq_samples/{sample}/{sample}.{size}_R1.fastq.gz',
        fastq_R2='data/rnaseq_samples/{sample}/{sample}.{size}_R2.fastq.gz'
    output:
        paired=temp('data/rnaseq_samples/{sample}/{sample}.{size}.paired.sam'),
        single=temp('data/rnaseq_samples/{sample}/{sample}.{size}.single.sam'),
    params: index='data/seq/dm6.full'
    threads: 8
    run:
        shell(
            'hisat2 '
            '-x {params.index} '
            '-1 {input.fastq_R1} '
            '-2 {input.fastq_R2} '
            '-p {threads} '
            '-S {output.paired}'
        )
        shell(
            'hisat2 '
            '-x {params.index} '
            '-U {input.fastq_R1} '
            '-p {threads} '
            '-S {output.single}'
        )

# ------------------------------------------------------------------------------
# HISAT2 outputs SAM but most tools use BAM
rule rnaseq_bam:
    input:
        paired=rules.hisat_align.output.paired,
        single=rules.hisat_align.output.single
    output:
        paired=temp('data/rnaseq_samples/{sample}/{sample}.{size}.paired.bam'),
        single=temp('data/rnaseq_samples/{sample}/{sample}.{size}.single.bam')
    run:
        shell('samtools view -Sb {input.paired} > {output.paired}')
        shell('samtools view -Sb {input.single} > {output.single}')

rule rnaseq_sortbam:
    input:
        paired=rules.rnaseq_bam.output.paired,
        single=rules.rnaseq_bam.output.single,
    output:
        paired='data/rnaseq_samples/{sample}/{sample}.{size}.paired.sorted.bam',
        single='data/rnaseq_samples/{sample}/{sample}.{size}.single.sorted.bam'
    shell:
        'samtools sort -m 4G {input.paired} > {output.paired} '
        '&& samtools sort {input.single} > {output.single} '


# ----------------------------------------------------------------------------
# Bowtie2 align.
#
rule bowtie2_align:
    input:
        index=expand('data/seq/dm6.full.{n}.bt2', n=range(1,2)),
        fastq_R1='data/chipseq_samples/{sample}/{sample}.{size}_R1.fastq.gz'
    output:
        single=temp('data/chipseq_samples/{sample}/{sample}.{size}.single.sam'),
    params: index='data/seq/dm6.full'
    threads: 8
    run:
        shell(
            'bowtie2 '
            '-x {params.index} '
            '-U {input.fastq_R1} '
            '-p {threads} '
            '-S {output.single}'
        )


# ------------------------------------------------------------------------------
# Bowtie2 outputs SAM but most tools use BAM
rule chipseq_bam:
    input:
        single=rules.bowtie2_align.output.single
    output:
        single=temporary('data/chipseq_samples/{sample}/{sample}.{size}.single.bam')
    run:
        shell('samtools view -Sb {input.single} > {output.single}')

rule chipseq_sortbam:
    input:
        single=rules.chipseq_bam.output.single,
    output:
        single='data/chipseq_samples/{sample}/{sample}.{size}.single.sorted.bam'
    shell:
        'samtools sort {input.single} > {output.single} '


# ----------------------------------------------------------------------------
# Bowtie2 align - cut-and-run
#
rule bowtie2_align_cutnrun:
    input:
        index=expand('data/seq/dm6.full.{n}.bt2', n=range(1,2)),
        fastq_R1='data/cutnrun_samples/{sample}/{sample}.{size}_R1.fastq.gz',
        fastq_R2='data/cutnrun_samples/{sample}/{sample}.{size}_R2.fastq.gz'
    output:
        paired=temporary('data/cutnrun_samples/{sample}/{sample}.{size}.paired.sam'),
    params: index='data/seq/dm6.full'
    threads: 8
    run:
        shell(
            'bowtie2 '
            '-x {params.index} '
            '-1 {input.fastq_R1} '
            '-2 {input.fastq_R2} '
            '-p {threads} '
            '-S {output.paired}'
        )


# ------------------------------------------------------------------------------
# Bowtie2 outputs SAM but most tools use BAM
rule cutnrun_bam:
    input:
        paired=rules.bowtie2_align_cutnrun.output.paired
    output:
        paired=temporary('data/cutnrun_samples/{sample}/{sample}.{size}.paired.bam')
    run:
        shell('samtools view -Sb {input.paired} > {output.paired}')

rule cutnrun_sortbam:
    input:
        paired=rules.cutnrun_bam.output.paired,
    output:
        paired='data/cutnrun_samples/{sample}/{sample}.{size}.paired.sorted.bam'
    shell:
        'samtools sort {input.paired} > {output.paired} '


# Previous iterations tried to subset the BAM by the limits and then convert to
# FASTQ. But since mates can span the limit boundaries, this can result in
# FASTQs with mismatched read counts.
#
# Still other iterations had tried to map to a restricted reference, but this
# still took a while, and any changes that would need to be made (e.g. in
# number of reads, or in the limits) would trigger a re-run from the beginning,
# including building the reference.
#
# Using seqtk here ensures that we're working at the read name level which,
# after all is what we care about for making a small FASTQ.
#
# Also note the parameters here for the maximum number of mapped and unmapped
# reads to include in the final small FASTQ.


rule rnaseq_small_fastq:
    input:
        bam='data/rnaseq_samples/{sample}/{sample}.full.paired.sorted.bam',
        full_fastq_R1=rules.download_rnaseq_fastqs.output.fastq_R1,
        full_fastq_R2=rules.download_rnaseq_fastqs.output.fastq_R2,
        limits=rules.limits.output
    output:
        mapped_names='data/rnaseq_samples/{sample}/{sample}.{size}.names.mapped.lst',
        unmapped_names='data/rnaseq_samples/{sample}/{sample}.{size}.names.unmapped.lst',
        R1='data/rnaseq_samples/{sample}/{sample}.{size}_R1.fastq',
        R2='data/rnaseq_samples/{sample}/{sample}.{size}_R2.fastq'
    run:
        mapped_n = mapped_n_config[wildcards.size]
        unmapped_n = unmapped_n_config[wildcards.size]

        shell(
            'samtools view -h -L {input.limits} {input.bam} '
            ' | samtools view -f 3 - '
            ' | cut -f1 '
            ' | sort -u '
            ' | head -n {mapped_n} > {output.mapped_names} '
        )

        shell(
            'samtools view -f 4 {input.bam} '
            ' | cut -f1 '
            ' | sort -u '
            ' | head -n {unmapped_n} > {output.unmapped_names} '
        )

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.mapped_names} '
            '> {output.R1}')

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.unmapped_names} '
            '>> {output.R1}')

        shell(
            'seqtk subseq {input.full_fastq_R2} {output.mapped_names} '
            '> {output.R2}')

        shell(
            'seqtk subseq {input.full_fastq_R2} {output.unmapped_names} '
            '>> {output.R2}')

rule chipseq_small_fastq:
    input:
        bam='data/chipseq_samples/{sample}/{sample}.full.single.sorted.bam',
        full_fastq_R1=rules.download_chipseq_fastqs.output.fastq_R1,
        limits=rules.limits.output
    output:
        uniquely_mapped_names='data/chipseq_samples/{sample}/{sample}.{size}.names.mapped.lst',
        multi_mapped_names='data/chipseq_samples/{sample}/{sample}.{size}.names.multi.lst',
        unmapped_names='data/chipseq_samples/{sample}/{sample}.{size}.names.unmapped.lst',
        R1='data/chipseq_samples/{sample}/{sample}.{size}_R1.fastq',
    run:
        mapped_n = mapped_n_config[wildcards.size]
        unmapped_n = unmapped_n_config[wildcards.size]
        multi_n = multimapped_n_config[wildcards.size]

        # For "unique" and "multimapping" we need mutually exclusive sets of
        # reads to avoid including reads twice, which will cause Picard tools
        # to fail. So for unique, we require MAPQ > 20 and no XS:i tag, and for
        # multimapping MAPQ < 20 plus an XS:i tag.
        #
        # Bowtie2 can report MAPQ > 20 but also include an XS:i tag, which
        # effectively means "this read maps well enough to call uniquely mapped
        # (MAPQ >20) but there's another read reported (with score XS:i). That
        # other read has a worse MAPQ that's poor enough such that we do not
        # consider this read a multimapper".
        shell(
            'samtools view -h -L {input.limits} {input.bam} '
            ' | samtools view -F 4 - '
            ''' | awk -F "\t" '{{if ($5 > 20) print $0}}' '''
            ' | grep -v "XS:i" '
            ' | cut -f1 '
            ' | sort -u '
            ' | head -n {mapped_n} > {output.uniquely_mapped_names} '
        )

        shell(
            'samtools view -h -L {input.limits} {input.bam} '
            ' | samtools view -F 4 - '
            ''' | awk -F "\t" '{{if ($5 < 20) print $0}}' '''
            ' | grep "XS:i" '
            ' | cut -f1 '
            ' | sort -u '
            ' | head -n {multi_n} > {output.multi_mapped_names} '
        )

        shell(
            ' samtools view -f 4 {input.bam} '
            ' | cut -f1 '
            ' | sort '
            ' | uniq '
            ' | head -n {unmapped_n} > {output.unmapped_names} '
        )

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.uniquely_mapped_names} '
            '> {output.R1} ')

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.multi_mapped_names} '
            '>> {output.R1} ')

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.unmapped_names} '
            '>> {output.R1} ')


rule cutnrun_small_fastq:
    input:
        bam='data/cutnrun_samples/{sample}/{sample}.full.paired.sorted.bam',
        full_fastq_R1=rules.download_cutnrun_fastqs.output.fastq_R1,
        full_fastq_R2=rules.download_cutnrun_fastqs.output.fastq_R2,
        limits=rules.limits.output
    output:
        uniquely_mapped_names='data/cutnrun_samples/{sample}/{sample}.{size}.names.mapped.lst',
        multi_mapped_names='data/cutnrun_samples/{sample}/{sample}.{size}.names.multi.lst',
        unmapped_names='data/cutnrun_samples/{sample}/{sample}.{size}.names.unmapped.lst',
        R1='data/cutnrun_samples/{sample}/{sample}.{size}_R1.fastq',
        R2='data/cutnrun_samples/{sample}/{sample}.{size}_R2.fastq',
    run:
        mapped_n = mapped_n_config[wildcards.size]
        unmapped_n = unmapped_n_config[wildcards.size]
        multi_n = multimapped_n_config[wildcards.size]

        # For "unique" and "multimapping" we need mutually exclusive sets of
        # reads to avoid including reads twice, which will cause Picard tools
        # to fail. So for unique, we require MAPQ > 20 and no XS:i tag, and for
        # multimapping MAPQ < 20 plus an XS:i tag.
        #
        # Bowtie2 can report MAPQ > 20 but also include an XS:i tag, which
        # effectively means "this read maps well enough to call uniquely mapped
        # (MAPQ >20) but there's another read reported (with score XS:i). That
        # other read has a worse MAPQ that's poor enough such that we do not
        # consider this read a multimapper".
        shell(
            'samtools view -h -L {input.limits} {input.bam} '
            ' | samtools view -F 4 - '
            ''' | awk -F "\t" '{{if ($5 > 20) print $0}}' '''
            ' | grep -v "XS:i" '
            ' | cut -f1 '
            ' | sort -u '
            ' | head -n {mapped_n} > {output.uniquely_mapped_names} '
        )

        shell(
            'samtools view -h -L {input.limits} {input.bam} '
            ' | samtools view -F 4 - '
            ''' | awk -F "\t" '{{if ($5 < 20) print $0}}' '''
            ' | grep "XS:i" '
            ' | cut -f1 '
            ' | sort -u '
            ' | head -n {multi_n} > {output.multi_mapped_names} '
        )

        shell(
            ' samtools view -f 4 {input.bam} '
            ' | cut -f1 '
            ' | sort '
            ' | uniq '
            ' | head -n {unmapped_n} > {output.unmapped_names} '
        )

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.uniquely_mapped_names} '
            '> {output.R1} ')
        shell(
            'seqtk subseq {input.full_fastq_R2} {output.uniquely_mapped_names} '
            '> {output.R2} ')

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.multi_mapped_names} '
            '>> {output.R1} ')
        shell(
            'seqtk subseq {input.full_fastq_R2} {output.multi_mapped_names} '
            '>> {output.R2} ')

        shell(
            'seqtk subseq {input.full_fastq_R1} {output.unmapped_names} '
            '>> {output.R1} ')
        shell(
            'seqtk subseq {input.full_fastq_R2} {output.unmapped_names} '
            '>> {output.R2} ')

rule gzipped_fastq:
    input:
        'data/{assay}_samples/{sample}/{sample}.{size}_R{N}.fastq'
    wildcard_constraints:
        size='small|tiny'
    output:
        'data/{assay}_samples/{sample}/{sample}.{size}_R{N}.fastq.gz'
    shell:
        'gzip {input}'

rule gzipped_gtf:
    input:
        rules.prep_small_gtf.input.gtf
    output: 'data/annotation/dm6.small.gtf.gz'
    shell:
        'gzip {input}'

# vim: ft=python
