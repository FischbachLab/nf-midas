# nf-midas

[![nf-midas](https://github.com/FischbachLab/nf-midas/actions/workflows/test.yaml/badge.svg)](https://github.com/FischbachLab/nf-midas/actions/workflows/test.yaml)

(fork of [midas_nextflow](https://github.com/kkerns85/midas_nextflow))

Running Midas using Nextflow
Optimized for Running on AWS Batch

> **Beta**
> This workflow is designed to help streamline bacterial metagenomic and metatranscriptomic data analysis using the Nextflow workflow manager. Trimmed, aligned, and filtered reads are passed on to MIDAS in order to assign taxonomy using the Phy-Eco single-copy marker gene set and pangenome database (v1.2) which provides species taxonomic level assignment, gene content, and strain level variation via single-nucleotide-polymorphism (SNP's) analysis for each metagenome.

This minimal test_data set was developed using Saccharibacteria Nanosynbacter lyticus HMT 952, formerly TM7x, a egnimatic member of the candidate phyla radiation (CPR), further emphazing the application of this workflow.

## Usage

### Test

```bash
aws batch submit-job \
    --job-name nf-midas-test-0302-2 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-midas,\
"--manifest","s3://genomics-workflow-core/Results/midas/00_TEST/00_seedfile/test.seedfile.csv"
```

### Actual sample (Paired-End)

```bash
aws batch submit-job \
    --job-name nf-midas-hComB6Gen \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-midas,\
"--project","hCom-B6Gen",\
"--manifest","s3://genomics-workflow-core/Results/midas/hCom-B6Gen/00_seedfile/hCom-B6Gen.seedfile.csv"
```

## Install Nextflow

- <https://www.nextflow.io/docs/latest/getstarted.html>

## Input

- Manifest file (S3://)
  - Path to paried ~~or single end~~ fastq or fastq.gz (R1 and R2) (S3://)
  - Metadata additional Columns with simple header format
- Work directory (S3://)
- Path to Midas Database (S3://)

## Output

By default, the outdir is `s3://genomics-workflow-core/Results/midas/<PROJECT>/`

- Species Analysis from MiDAS
- Gene Analysis from MiDAS
- SNP Analysis from MiDAS
- Merged files for Species, Genes, and SNP analysis

## Options

- --outdir     Folder to place analysis outputs (default ./midas)
- --species_cov       Coverage (depth) threshold for species inclusion (default: 3.0)
- --merge_sample_depth  Corresponds to the --sample_depth parameter in the merge_midas.py command (default: 1.0)

## `--manifest`

- The manifest is a CSV with a header indicating which samples correspond to which files.
- The file must contain a column `sampleName`. This value can not be repeated.
- Reads are specified by columns, `R1` and `R2`.

## Software Incorporated

- MiDAS Docker: quay.io/biocontainers/midas:1.3.2--pyh5e36f6f_6
- [MiDAS](https://github.com/snayfach/MIDAS)
- [Nextflow](https://www.nextflow.io)

## Developers

- [Sunit Jain](www.sunitjain.com)
- Xiandong Meng
