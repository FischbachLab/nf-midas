
process midas {
    tag "$sampleName"
    container params.docker_container_midas
    label "mem_veryhigh"
    publishDir "${params.outdir}/${params.project}/${sampleName}"

    input:
    tuple val(sampleName), path("${sampleName}.R1.fastq.gz"), path ("${sampleName}.R2.fastq.gz")

    output:
    path "${sampleName}.species.tar.gz", emit: species_ch
    path "${sampleName}.genes.tar.gz", optional: true, emit: gene_ch
    path "${sampleName}.snps.tar.gz", optional: true, emit: snps_ch
    path "job.*", optional: true
    path "*.log.txt", optional: true

    """
    ls -lhtra
    export FWD_FASTQ="${sampleName}.R1.fastq.gz"; \
    export REV_FASTQ="${sampleName}.R2.fastq.gz"; \
    export DB=${params.db_midas}; \
    export CORES=${task.cpus}; \
    export SAMPLE_NAME=${sampleName}; \
    run_midas.sh
    """
}

process midas_merge_species {
    container params.docker_container_midas
    label "mem_veryhigh"

    publishDir "${params.outdir}/${params.project}/merged_species/"

    input:
    path species_tar_list
    // file DB from file(params.db_midas)

    output:
    path "SPECIES/*"

"""
#!/bin/bash -x
set -e
ls -lahtr
# Keep track of the folders created while unpacking input files
input_string=""
echo "Unpacking all of the input files"
for tarfile in ${species_tar_list}; do
    echo "Making sure that \$tarfile was downloaded correctly"
    [[ -s \$tarfile ]]
    echo "Unpacking \$tarfile"
    tar xzvf \$tarfile
    # Add this folder to the input string
    input_string="\$input_string,\$( echo \$tarfile | sed 's/.species.tar.gz//' )"
    echo "Updated input string: \$input_string"
done
# Remove the leading comma from the input string
input_string=\$( echo \$input_string | sed 's/^,//' )
if [ -z \$input_string ]; then
    echo "[SPECIES]: Input string is empty. Nothing to do here."
    mkdir -p SPECIES
    touch SPECIES/NOT_ENOUGH_DATA
else
    echo "Merging species results"
    merge_midas.py \
        species \
        SPECIES \
        -i \$input_string \
        -t list \
        -d ${params.db_midas} \
        --sample_depth ${params.merge_sample_depth}
    touch SPECIES/DONE
fi
echo "Done merging data"
ls -lahtr SPECIES
echo "Compressing output files"
find SPECIES -type f | xargs gzip
echo "Done"
"""
}


process midas_merge_genes {
    container params.docker_container_midas
    label "mem_veryhigh"
    publishDir "${params.outdir}/${params.project}/merged_genes/"

    input:
    path genes_tar_list

    output:
    path "GENES/*"

"""
#!/bin/bash -x
set -e
ls -lahtr
# Keep track of the folders created while unpacking input files
input_string=""
echo "Unpacking all of the input files"
for tarfile in ${genes_tar_list}; do
    echo "Making sure that \$tarfile was downloaded correctly"
    [[ -s \$tarfile ]]
    echo "Unpacking \$tarfile"
    tar xzvf \$tarfile
    # Add this folder to the input string
    input_string="\$input_string,\$( echo \$tarfile | sed 's/.genes.tar.gz//' )"
    echo "Updated input string: \$input_string"
done
# Remove the leading comma from the input string
input_string=\$( echo \$input_string | sed 's/^,//' )
if [ -z \$input_string ]; then
    echo "[GENES]: Input string is empty. Nothing to do here."
    mkdir -p GENES
    touch GENES/NOT_ENOUGH_DATA
else
    echo "Merging gene results"
    merge_midas.py \
        genes \
        GENES \
        -i \$input_string \
        -t list \
        -d ${params.db_midas} \
        --sample_depth ${params.merge_sample_depth}
    touch GENES/DONE
fi
echo "Done merging data"
ls -lahtr GENES
echo "Compressing output files"
find GENES -type f | xargs gzip
echo "Done"
"""
}

process midas_merge_snps {
    container params.docker_container_midas
    label "mem_veryhigh"
    publishDir "${params.outdir}/${params.project}/merged_snps/"

    input:
    path snps_tar_list
    // path DB from file(params.db_midas)

    output:
    path "SNPS/*"

"""
#!/bin/bash -x
set -e
ls -lahtr
# Keep track of the folders created while unpacking input files
input_string=""
echo "Unpacking all of the input files"
for tarfile in ${snps_tar_list}; do
    echo "Making sure that \$tarfile was downloaded correctly"
    [[ -s \$tarfile ]]
    echo "Unpacking \$tarfile"
    tar xzvf \$tarfile
    # Add this folder to the input string
    input_string="\$input_string,\$( echo \$tarfile | sed 's/.snps.tar.gz//' )"
    echo "Updated input string: \$input_string"
done
# Remove the leading comma from the input string
input_string=\$( echo \$input_string | sed 's/^,//' )

if [ -z \$input_string ]; then
    echo "[SNPS]: Input string is empty. Nothing to do here."
    mkdir -p SNPS
    touch SNPS/NOT_ENOUGH_DATA
else
    echo "Merging snps results"
    merge_midas.py \
        snps \
        SNPS \
        -i \$input_string \
        -t list \
        -d ${params.db_midas} \
        --sample_depth ${params.merge_sample_depth}
    touch SNPS/DONE
fi
echo "Done merging data"
ls -lahtr SNPS
echo "Compressing output files"
find SNPS -type f | xargs gzip
echo "Done"
"""
}
