#!/bin/bash

set -euo pipefail

ACTUAL_BINARY_LOC="/usr/local/bin"
EXPECTED_BINARY_LOC="/usr/local/lib/python3.9/site-packages/bin/Linux"

mkdir -p ${EXPECTED_BINARY_LOC}

# From Midas, utility.py (https://github.com/snayfach/MIDAS/blob/7ac8d6e6515b5b780132494102e4d1883fbf407b/midas/utility.py#L114)
# args['hs-blastn'] = '/'.join([main_dir, 'bin', platform.system(), 'hs-blastn'])
# args['bowtie2-build'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2-build'])
# args['bowtie2'] = '/'.join([main_dir, 'bin', platform.system(), 'bowtie2'])
# args['samtools'] = '/'.join([main_dir, 'bin', platform.system(), 'samtools'])
ln -s ${ACTUAL_BINARY_LOC}/hs-blastn ${EXPECTED_BINARY_LOC}/hs-blastn
ln -s ${ACTUAL_BINARY_LOC}/bowtie2* ${EXPECTED_BINARY_LOC}/
ln -s ${ACTUAL_BINARY_LOC}/samtools ${EXPECTED_BINARY_LOC}/samtools