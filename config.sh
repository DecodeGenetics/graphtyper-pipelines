#!/usr/bin/env bash

set -e

# Files #
## (Required) Path to your indexed reference genome, e.g. /nfs/data/GRCh38.fa
GENOME=

## (Optional) Location of indexed bgzipped VCF file to use to initialize the graph structure. Should not include SVs (they go below)
VCF=

## (Optional) VCF with candicate SVs
SV_VCF=

# Program binaries #
## (Required) graphtyper binary, e.g. /usr/bin/graphtyper
GRAPHTYPER=$(type -P graphtyper || true)

## (Optional in SNP/indel calling, required in SV calling) bamShrink binary
BAMSHRINK=$(type -P bamShrink || true)

## (Required) tabix binary, e.g. /usr/bin/tabix
TABIX=$(type -P tabix || true)

## (Required) samtools binary, e.g. /usr/bin/samtools
SAMTOOLS=$(type -P samtools || true)

## (Required) GNU parallel binary, e.g. /usr/bin/parallel
PARALLEL=$(type -P parallel || true)


# Directories #
## Format of the temporary directory. Change this if you cannot use /tmp/
TMP_FORMAT="/tmp/graphtyper_calling.XXXXXX"

## Top directory, you should probably not change this unless you know what you are doing
TOP_DIR="$(realpath $(dirname ${BASH_SOURCE[0]}))"


# Region/job parameters #
## (optional) If you want only to genotype a specific region, you can define its start position here.
## E.g. "chr1:6180000". Leave empty to call entire genome.
REGION_START=

## (optional) If you want to genotype a list of regions you can specify a file with them here.
## The regions should be one per line and in the format: chr1:
REGION_FILE=

## Size of the region each job will cover. Should probably never be more than 1 MB
REGION_SIZE=1000000 # 1 MB

## Number of threads each job will be allocated
NUM_THREADS=1

## Number of threads in SV calling
NUM_THREADS_SV_CALLING=24

## Number of slices to run at the same time in each region/job.
NUM_SLICES_RUNNING=1

## Whether or not the temporary directories should be cleaned after genotyping
CLEAN_UP=1

## Chromosomes to call
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 \
chr19 chr20 chr21 chr22"

## Set to 1 if there should only one genotyping iterations (only activated with a VCF file is given)
GENOTYPE_ONLY=0

# Call parameters #
## Number of threads will use to genotype each slice
## 1 is recommended for most cases, but you can use more if you want to reduce the total time of each slice.
GRAPHTYPER_THREADS=1

## Graphtyper call options for large sample sizes
GRAPHTYPER_COMMON_OPTS="--minimum_variant_support=5 --minimum_variant_support_ratio=0.35"

## Graphtyper maximum allowed number of files open at the same time
MAX_FILES_OPEN=9800

## Number of bases in each slice.
SLICE_SIZE=50000

## Number of bases padded around slices
PAD_SIZE=200

## Extra padding applied in SV calling
EXTRA_SV_PADDING=200000

## Results of SV genotyping
SV_RESULTS_DIR=sv_results

## Set to 1 if the first sample should be appended to the SV results directory name
SV_RESULTS_DIR_SUFFIX_SAMPLE1=0

# Read my_config.sh if available (note: that file should not be version controlled) #
if [[ -f my_config.sh ]]
then
  source my_config.sh
fi


# Check for problems#
if [[ ! -f $GRAPHTYPER ]]; then echo "Graphtyper was not found (GRAPHTYPER was to '$GRAPHTYPER')."; exit 1; fi
if [[ ! -f $TABIX ]]; then echo "Tabix was not found (TABIX was set to '$TABIX')."; exit 1; fi
if [[ ! -f $SAMTOOLS ]]; then echo "SAMtools was not found (SAMTOOLS was set to '$SAMTOOLS')."; exit 1; fi
if [[ ! -f $PARALLEL ]]; then echo "GNU parallel was not found (PARALLEL was set to '$PARALLEL')."; exit 1; fi
if [[ ! -f $GENOME ]]; then echo "Reference genome was not found (GENOME was set to '$GENOME')."; exit 1; fi
if [[ ! -z $VCF ]] && [[ ! -f $VCF ]]; then echo "VCF was not found (VCF was set to '$VCF')."; exit 1; fi

# if the region size is not dividable by the slice size we need to make it larger
rem=$((REGION_SIZE % SLICE_SIZE))

if [[ $rem -ne 0 ]]; then
  REGION_SIZE=$((REGION_SIZE - rem + SLICE_SIZE))
fi
