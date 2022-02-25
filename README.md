# nf-midas

(fork of [midas_nextflow](https://github.com/kkerns85/midas_nextflow))

Running Midas using Nextflow
Optimized for Running on AWS Batch

> **Beta**
> This workflow is designed to help streamline bacterial metagenomic and metatranscriptomic data analysis using the Nextflow workflow manager. Raw reads (fastq or fastq.gz) are initailly passed through Kneaddata which initially trims raw reads via Trimmomatic (v. 0.33) and then aligns them via Bowtie2 (v. >= 2.2) and finally removes bacterial ribosomal RNA using the Silva rRNA database (v. 128). Trimmed, aligned, and filtered reads are then passed on to MIDAS in order to assign taxonomy using the Phy-Eco single-copy marker gene set and pangenome database (v1.2) which provides species taxonomic level assignment, gene content, and strain level variation via single-nucleotide-polymorphism (SNP's) analysis for each metagenome.

This minimal test_data set was developed using Saccharibacteria Nanosynbacter lyticus HMT 952, formerly TM7x, a egnimatic member of the candidate phyla radiation (CPR), further emphazing the application of this workflow.

## Install Nextflow

- <https://www.nextflow.io/docs/latest/getstarted.html>

## [OPTIONAL] Download/Update the Kneaddata database

```{bash}
mkdir -p /mnt/efs/databases/Biobakery/kneaddata
cd /mnt/efs/databases/Biobakery/kneaddata
docker container run \
  --rm \
  --volume $PWD:$PWD \
  --workdir $PWD \
  biobakery/workflows:3.0.0.a.7 \
    kneaddata_database \
      --download human_genome \
      bowtie2 .
```

As of `2022-02-23` here are the databases that were available:

```{bash}
KneadData Databases ( database : build = location )
human_genome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
human_genome : bmtagger = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_BMTagger_v0.1.tar.gz
human_transcriptome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz
ribosomal_RNA : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2.tar.gz
mouse_C57BL : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz
```

## Input

- Manifest file (Local)
  - Path to paried or single end fastq or fastq.gz (R1 and R2) (Local or S3://)
  - Metadata additional Columns with simple header format
- Work directory (Local or S3://)
- Path to Midas Database (Local or S3://)

## Output

- Species Analysis from MiDAS
- Gene Analysis from MiDAS
- SNP Analysis from MiDAS
- Merged files for Species, Genes, and SNP analysis

## Options

- --single            Run single end reads
- --profile           local or AWS batch
- --output_folder     Folder to place analysis outputs (default ./midas)
- --output_prefix     Text used as a prefix for output files (default: midas)
- --species_cov       Coverage (depth) threshold for species inclusion (default: 3.0)
- --merge_sample_depth  Corresponds to the --sample_depth parameter in the merge_midas.py command (default: 1.0)

## Manifest

- The manifest is a CSV with a header indicating which samples correspond to which files.
- The file must contain a column `sampleName`. This value can not be repeated.
- Reads are specified by columns, `R1` and `R2`.

(Fix these; Seems contradictory)

- Data is only accepted as paired reads.
- If you specify --single, then only data from `R1` will be used

## Software Incorporated

More information is available at:

Kneaddata: <http://huttenhower.sph.harvard.edu/kneaddata>

Trimmomatic: <http://www.usadellab.org/cms/?page=trimmomatic>

Bowtie2: <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>

MiDAS: <https://github.com/snayfach/MIDAS>

MiDAS Docker: <https://github.com/FredHutch/docker-midas>

Nextflow: <https://www.nextflow.io>

## Developers

- [Sunit Jain](www.sunitjain.com)
- Xiandong Meng
