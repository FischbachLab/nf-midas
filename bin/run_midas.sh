#!/bin/bash

set -euo pipefail

# FWD_FASTQ
# DB
CORES=${CORES:1}
REV_FASTQ=${REV_FASTQ:''}

exit_handler () {
    # accepts 2 variables. 1) Exit status 2) Message
    local EXIT_STATUS=$1
    local MESSAGE="${2}"
    DURATION=$((SECONDS - START_TIME))
    hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
    local STATUS_FILE="job.complete"
    if [ "${EXIT_STATUS}" -ne "0" ]; then
        STATUS_FILE="job.failed"
    fi

    echo "${MESSAGE}" > "${STATUS_FILE}"
    printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs >> "${STATUS_FILE}"
    echo "Live long and prosper" >> "${STATUS_FILE}"
    ############################ PEACE! ################################
    # exit ${EXIT_STATUS}
}

MIDAS_SPECIES_PARAMS="-t ${CORES} -d ${DB} --remove_temp"
### Single end
if [[ -n ${FWD_FASTQ} ]]; then
    MIDAS_SPECIES_PARAMS="${MIDAS_SPECIES_PARAMS} -1 ${FWD_FASTQ}"
else
    # Probably 
    exit_handler 1 "[FATAL] Missing FWD fastq file."
fi

### Paired
if [[ -n ${REV_FASTQ:-} ]]; then
    MIDAS_SPECIES_PARAMS="${MIDAS_SPECIES_PARAMS} -2 ${REV_FASTQ}"
else
    # Probably 
    exit_handler 1 "[FATAL] Missing REV fastq file."
fi

## Set up parameters for Midas Genes and SNP
MIDAS_GENES_PARAMS=${MIDAS_SPECIES_PARAMS}
MIDAS_SNP_PARAMS=${MIDAS_SPECIES_PARAMS}

## Execution commands for MIDAS jobs
MIDAS_SPECIES="run_midas.py species ${SAMPLE_NAME} ${MIDAS_SPECIES_PARAMS} | tee -a midas_species.log.txt"
MIDAS_GENES="run_midas.py genes ${SAMPLE_NAME} ${MIDAS_GENES_PARAMS} | tee -a midas_genes.log.txt"
MIDAS_SNP="run_midas.py snps ${SAMPLE_NAME} ${MIDAS_SNP_PARAMS} | tee -a midas_snp.log.txt"

####################### MIDAS - Species ############################

if eval "${MIDAS_SPECIES}"; then
    echo "[$(date)] MIDAS Species complete." | tee -a midas_species.log.txt
    echo "Tarring up species results"
    tar cvzf ${SAMPLE_NAME}.species.tar.gz ${SAMPLE_NAME}/species/*
else
    exit_handler 1 "[FATAL] MIDAS Species failed."
fi

######################## MIDAS - Genes #############################

if eval "${MIDAS_GENES}"; then
    echo "[$(date)] MIDAS Genes complete." | tee -a midas_genes.log.txt
    echo "Tarring up gene results"
    tar cvzf ${SAMPLE_NAME}.genes.tar.gz ${SAMPLE_NAME}/genes/*
else
    exit_handler 0 "[FATAL] MIDAS Genes failed."
fi

######################## MIDAS - SNPs #############################

if eval "${MIDAS_SNP}"; then
    echo "[$(date)] MIDAS SNP complete." | tee -a midas_snps.log.txt
    echo "Tarring up SNP results"
    tar cvzf ${SAMPLE_NAME}.snps.tar.gz ${SAMPLE_NAME}/snps/*
echo "Done"
else
    exit_handler 0 "[FATAL] MIDAS SNPs failed."
fi

######################### HOUSEKEEPING #############################
exit_handler 0 "Success!"