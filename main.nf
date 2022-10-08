#!/usr/bin/env nextflow
nextflow.enable.dsl=1
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Analyze microbial pan-genomes using MIDAS
    Usage:
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)
      --db_midas            Folder containing the MIDAS database
      --project             Broad project/study this data belongs to.
      --prefix              Specific version of the data being analyzed. Todays date in YYYYMMDD format is acceptable if prefix is unknown.
    Options:
      --outdir       Folder to place analysis outputs (default ./midas)
      --species_cov         Coverage (depth) threshold for species inclusion (default: 3.0)
      --merge_sample_depth  Corresponds to the --sample_depth parameter in the merge_midas.py command (default: 1.0)
    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `sampleName`. This value can not be repeated.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If you specify --single, then only data from `R1` will be used
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime, or omits the --manifest
if (params.help || params.manifest == null || params.db_midas == null || params.db_knead == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make sure that the manifest file can be found
if (file(params.manifest).isEmpty()){

    // Print a helpful log message
    log.info"""
    Cannot find the file specified by --manifest ${params.manifest}
    """.stripIndent()

    // Exit out and do not run anything else
    exit 0
}

// Make sure that the Midas database file can be found
// if (file(params.db_midas).isEmpty()){

//     // Print a helpful log message
//     log.info"""
//     Cannot find the file specified by --db_midas ${params.db_midas}
//     """.stripIndent()

//     // Exit out and do not run anything else
//     exit 0
// }

// Make sure that the Kneaddata database file can be found
// if (file(params.db_knead).isEmpty()){

//     // Print a helpful log message
//     log.info"""
//     Cannot find the file specified by --db_knead ${params.db_knead}
//     """.stripIndent()

//     // Exit out and do not run anything else
//     exit 0
// }

// Channel.fromPath( "${params.db_knead}", checkIfExists: true ).set { db_knead_ch }

// Parse the manifest CSV
// Along the way, make sure that the appropriate columns were provided
if (params.single){
    fastq_ch = Channel.from(
        file(
            params.manifest
        ).splitCsv(
            header: true,
            sep: ","
        )
    ).filter {
        r -> (r.sampleName != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'sampleName' column: ${params.manifest}"
    }.filter {
        r -> (r.R1 != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'R1' column: ${params.manifest}"
    }.map {
        r -> [r["sampleName"], [file(r["R1"])]]
    }
} else {
    fastq_ch = Channel.from(
        file(
            params.manifest
        ).splitCsv(
            header: true,
            sep: ","
        )
    ).filter {
        r -> (r.sampleName != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'sampleName' column: ${params.manifest}"
    }.filter {
        r -> (r.R1 != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'R1' column: ${params.manifest}"
    }.filter {
        r -> (r.R2 != null)
    }.ifEmpty {
        exit 1, "Cannot find values in the 'R2' column: ${params.manifest}"
    }.map {
        r -> [r["sampleName"], [file(r["R1"]), file(r["R2"])]]
    }
}


//Creates working dir
workingpath = params.outdir + "/" + params.project
workingdir = file(workingpath)

if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() )     {
        exit 1, "Cannot create working directory: $workingpath"
    } 
}    

if(params.prefix){
    workingpath = workingpath + "/" + params.prefix
}

process midas {
    tag "$sampleName"
    container params.docker_container_midas
    label "mem_veryhigh"
    publishDir "${workingpath}/${sampleName}"

    input:
    tuple val(sampleName), file("${sampleName}.R*.fastq.gz") from fastq_ch

    output:
    file "${sampleName}.species.tar.gz" into species_ch
    file "${sampleName}.genes.tar.gz" optional true into gene_ch
    file "${sampleName}.snps.tar.gz" optional true into snps_ch
    file "job.*" optional true
    file "*.log.txt" optional true

    """
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
    publishDir "${workingpath}"

    input:
    file species_tar_list from species_ch.toSortedList()
    // file DB from file(params.db_midas)

    output:
    file "SPECIES/*"

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
if [ -z $input_string ]; then
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
    publishDir "${workingpath}"

    input:
    file genes_tar_list from gene_ch.toSortedList()

    output:
    file "GENES/*"

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
if [ -z $input_string ]; then
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
    publishDir "${workingpath}"

    input:
    file snps_tar_list from snps_ch.toSortedList()
    // file DB from file(params.db_midas)

    output:
    file "SNPS/*"

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

if [ -z $input_string ]; then
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
