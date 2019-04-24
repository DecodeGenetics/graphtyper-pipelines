#!/usr/bin/env bash

set -e
set -o pipefail


bam=$1
CONFIG=$2

if [[ -z "$bam" ]]
then
  echo "Usage: $0 <BAM/bamlist> [CONFIG]"
  exit 1
fi

TOP_DIR="$(realpath $(dirname ${BASH_SOURCE[0]}))"

if [[ -z "$CONFIG" ]]
then
  CONFIG=config.sh
fi

source $CONFIG

# Concatenate haplotypes
for chrom in $CHROMOSOMES
do
  # Skip check if there is already a VCF with haplotypes
  if [[ -f haps/${chrom}.vcf.gz ]] && [[ -f haps/${chrom}.vcf.gz.tbi ]]; then continue; fi

  start=1
  end=$((start - 1 + REGION_SIZE))
  CHROM_SIZE=`grep -w "^${chrom}" $GENOME.fai | cut -f2`
  end=$((end>CHROM_SIZE?CHROM_SIZE:end)) # Never expand end further than the total length of the chromosome

  while [[ $start -lt $CHROM_SIZE ]]
  do
    region_id=`printf "%09d" $start`"-"`printf "%09d" $end`

    # Only add this region if the output file is missing
    if [[ ! -f "haps/${chrom}/${region_id}.vcf.gz" ]]; then
      echo "ERROR: Missing input file 'haps/${chrom}/${region_id}.vcf.gz'. Aborting." >&2
      exit 1
      #echo "set -e; set -o pipefail; ./node_script.sh $CONFIG $bam ${chrom}:${start}"
    fi

    echo "haps/${chrom}/${region_id}.vcf.gz"
    start=$((end + 1))
    end=$((start - 1 + REGION_SIZE))
    end=$((end>CHROM_SIZE?CHROM_SIZE:end)) # Never expand end further than the total length of the chromosome
  done > /tmp/${chrom}_filelist

  $VT cat -s -L /tmp/${chrom}_filelist -o haps/${chrom}.vcf.gz
  $TABIX haps/${chrom}.vcf.gz
done

# echo commands to run on a cluster
SV_REGION_SIZE=1000000 # 1MB

for chrom in $CHROMOSOMES
do
  start=1
  end=$((start - 1 + SV_REGION_SIZE))
  CHROM_SIZE=`grep -w "^${chrom}" $GENOME.fai | cut -f2`
  end=$((end>CHROM_SIZE?CHROM_SIZE:end)) # Never expand end further than the total length of the chromosome

  while [[ $start -lt $CHROM_SIZE ]]
  do
    region_id=`printf "%09d" $start`"-"`printf "%09d" $end`

    # Only add this region if the output file is missing
    if [[ ! -f "sv_results/${chrom}/${region_id}.vcf.gz" ]]; then
      echo "set -e; set -o pipefail; ./sv_node_script.sh $CONFIG $bam ${chrom}:${start}"
    fi

    start=$((end + 1))
    end=$((start - 1 + SV_REGION_SIZE))
    end=$((end>CHROM_SIZE?CHROM_SIZE:end)) # Never expand end further than the total length of the chromosome
  done
done
