#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { midas; midas_merge_species;} from './modules/midas'
include { midas_merge_genes; midas_merge_snps;} from './modules/midas'
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
      --outdir              Folder to place analysis outputs (default ./midas)
      --species_cov         Coverage (depth) threshold for species inclusion (default: 3.0)
      --merge_sample_depth  Corresponds to the --sample_depth parameter in the merge_midas.py command (default: 1.0)
      --genes               if enabling the merged gene analysis (default: false)
      --snps                if enabling the merged snp analysis (default: false)
    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `sampleName`. This value can not be repeated.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If you specify --single, then only data from `R1` will be used
    """.stripIndent()
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
 
 

workflow {
   
// Show help message if the user specifies the --help flag at runtime, or omits the --manifest
if (params.help || params.manifest == null || params.db_midas == null || params.db_knead == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Parse the manifest CSV
// Along the way, make sure that the appropriate columns were provided
/*
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

*/

    fastq_ch = Channel
            .fromPath(params.manifest)
            .ifEmpty { exit 1, "Cannot find any manifest file matching: ${params.manifest}" }        
            .splitCsv(header: ['sampleName', 'R1', 'R2'], sep: ',', skip: 1)
            .map{ row -> tuple(row.sampleName, row.R1, row.R2)}

    fastq_ch | midas

    midas_merge_species( midas.out.species_ch.toSortedList() )

    if (params.genes){
        midas_merge_genes( midas.out.genes_ch.toSortedList() )
    }

    if (params.snps){
        midas_merge_snps( midas.out.snps_ch.toSortedList() )
    }

}