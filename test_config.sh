#!/usr/bin/bash
set -e
set -o pipefail

# Files #
## (Required) Path to your indexed reference genome, e.g. /nfs/data/GRCh38.fa
GENOME=test/reference.fa

## (Optional) Location of indexed bgzipped VCF file to use to initialize the graph structure.
## Only needed if "INITIALIZE_GRAPH_WITH_VCF" is not 0. E.g. /nfs/data/dbSNP_common.vcf.gz
VCF=


# Program binaries #
## graphtyper binary, e.g. /usr/bin/graphtyper
GRAPHTYPER=$(type -P graphtyper)

## vt binary (https://github.com/atks/vt), e.g. /usr/bin/vt
VT=$(type -P vt)

## tabix binary, e.g. /usr/bin/tabix
TABIX=$(type -P tabix)

## samtools binary, e.g. /usr/bin/samtools
SAMTOOLS=$(type -P samtools)

## GNU parallel binary, e.g. /usr/bin/parallel
PARALLEL=$(type -P parallel)


# Directories #
## Format of the temporary directory. Change this if you cannot use /tmp/
TMP_FORMAT="/tmp/graphtyper_calling.XXXXXX"

## Top directory, you should probably not change this unless you know what you are doing
TOP_DIR="$(realpath $(dirname ${BASH_SOURCE[0]}))"

## Where the final results should go
RESULTS="${TOP_DIR}/test/results"


# Region/job parameters #
## (optional) If you want only to genotype a specific region, you can define its start position here.
## E.g. "chr1:6180000". Leave empty to call entire genome.
REGION_START=20:500

## Size of the region each job will cover
REGION_SIZE=9000

## Number of threads each job will be allocated
NUM_THREADS=4

## Number of slices to run at the same time in each region/job.
## Because of I/O operations, you may want to use a value which is more than your total thread count to fully utilize your CPU power.
NUM_SLICES_RUNNING=8

## Whether or not the temporary directories should be cleaned after genotyping
CLEAN_UP=1

## Chromosomes to call
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20\
 chr21 chr22 chrX"


# Call parameters #
## Should Graphtyper initialize the graph using a VCF file (specified with the "VCF" variable)
INITIALIZE_GRAPH_WITH_VCF=0

## Number of threads will use to genotype each slice
## 1 is recommended for most cases, but you can use more if you want to reduce the total time of each slice.
GRAPHTYPER_THREADS=1

## Number of samples to consider it to be a small sample size
SMALL_SAMPLE_SIZE=2

## Graphtyper call options for small sample sizes
GRAPHTYPER_SMALL_SAMPLE_SIZE_OPTS="--threads=${GRAPHTYPER_THREADS} --minimum_variant_support=3 --minimum_variant_support_ratio=0.15"

## Graphtyper call options for large sample sizes
GRAPHTYPER_POPULATION_OPTS="--threads=${GRAPHTYPER_THREADS} --minimum_variant_support=4 --minimum_variant_support_ratio=0.18"

## Number of bases in each slice.
SLICE_SIZE=4500

## Number of bases padded around slicess
PAD_SIZE=200


# Read my_config.sh if available (note: that file should not be version controlled) #
if [[ -f my_config.sh ]]
then
  source my_config.sh
fi


# Check for problems#
if [[ ! -f $GRAPHTYPER ]]; then echo "Graphtyper was not found (GRAPHTYPER was to '$GRAPHTYPER')."; exit 1; fi
if [[ ! -f $VT ]]; then echo "VT was not found (VT was set to '$VT')."; exit 1; fi
if [[ ! -f $TABIX ]]; then echo "Tabix was not found (TABIX was set to '$TABIX')."; exit 1; fi
if [[ ! -f $SAMTOOLS ]]; then echo "SAMtools was not found (SAMTOOLS was set to '$SAMTOOLS')."; exit 1; fi
if [[ ! -f $PARALLEL ]]; then echo "GNU parallel was not found (PARALLEL was set to '$PARALLEL')."; exit 1; fi
if [[ ! -f $GENOME ]]; then echo "Reference genome was not found (GENOME was set to '$GENOME')."; exit 1; fi
if [[ ! -f $VCF ]] && [[ $INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then echo "VCF was not found (VCF was set to '$VCF')."; exit 1; fi

# if the region size is not dividable by the slice size we need to make it larger
rem=$((REGION_SIZE % SLICE_SIZE))

if [[ $rem -ne 0 ]]; then
  REGION_SIZE=$((REGION_SIZE - rem + SLICE_SIZE))
fi
