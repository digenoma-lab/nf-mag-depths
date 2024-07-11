![Language: Nextflow](https://img.shields.io/badge/Language-Nextflow-green.svg)

# nf-mag-depths

A [nextflow](https://www.nextflow.io) pipeline to calculate depth of coverage from a metagenomic set of bins.

## Usage

First, we need our initial set of bins to be pre-processed and merged into a single FASTA file (```.fna```, ```.fa```). Also, a comma-separated input file declaring information about samples and their raw short-reads (```.fastq``` or ```.fq``` -compressed- files).

```
#usage/help:
nextflow run nf-mag-depths/main.nf -c nf-mag-depths/mag_depths.config --help
```

The pipeline uses [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) to index the reference assembly and to map it back to the raw reads for each metagenomic sample. After, it uses samtools to generate a sorted ```.bam``` file, which will be the input for ```jgi_summarize_bam_contig_depths``` (a tool from the [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) MAG binner suite) to calculate the coverage or depths of all bins from the reference across all samples.

## Input file

See the ```input.csv``` file as an example of how to declare the parameters regarding sample names and short-reads paths.

```
sample,read1,read2
S1,reads/S1.R1.fq.gz,reads/S1.R2.fq.gz
S2,reads/S2.R1.fq.gz,reads/S2.R2.fq.gz
S3,reads/S3.R1.fq.gz,reads/S3.R2.fq.gz
S7,reads/S7.R1.fq.gz,reads/S7.R2.fq.gz
S15,reads/S15.R1.fq.gz,reads/S15.R2.fq.gz
S17,reads/S17.R1.fq.gz,reads/S17.R2.fq.gz
```

## Conda environment

The needed environment for running the pipeline is declared in the ```env.yml``` file (included).

```
micromamba create -f env.yml
micromamba activate mag_depths
micromamba deactivate
```

## Config file

Last, a config file for nextflow is necessary for the execution, in which we declare global variables along with processes parameters. See the file: ```mag_depths.config```.

## Run

Once we have set all execution parameters, files, paths and other configurations, we can run the pipeline as follows:

```
nextflow -bg run nf-mag-depths/main.nf -c nf-mag-depths/mag_depths.config --input input.csv -profile singularity {-resume}
```

## Software versions

| Software | Version |
| :------- | :------ |
| Nextflow | 23.04.2 |
| bwa-mem2 | 2.2.1   |
| samtools | 1.16.1  |
| MetaBAT2 | 2.15    |
