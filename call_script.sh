#!/usr/bin/env bash

set -e
set -o pipefail

##
# help functions
##
get_hash()
{
  # argument $1 is a bgzipped vcf file
  hash=$(zgrep -v "^#" $1 | cut -f3-5 | md5sum | cut -d' ' -f1 || true)
  echo $hash
}


##
# Start of script
##

CONFIG=$1
bamlist=$2
region=$3 # E.g. chr21:1
TMP=`dirname $bamlist`

source $CONFIG

chrom=`echo $region | cut -d":" -f1`
b=`echo $region | cut -d":" -f2`
region_id=`printf "%09d" $b`"-"`printf "%09d" $((b - 1 + SLICE_SIZE))`

TMPR="$TMP/${chrom}/${region_id}"

mkdir --parents $TMPR/bams $TMPR/it{1,2,3,4,5}

# Clean up after calling
if [[ $CLEAN_UP -ne 0 ]]; then
  trap "rm -r -f $TMPR; exit 1" 1 2 15
fi

NUM_SAMPLES=`cat $bamlist | wc -l`
test `tail -c1 $bamlist` && NUM_SAMPLES=$((NUM_SAMPLES + 1))

# Graphtyper settings
GRAPH=$TMPR/graph
GT_LOG=$TMPR/gt_log

UNPADDED_REGION="${chrom}:${b}-$((b - 1 + SLICE_SIZE))"
PADDED_REGION="${chrom}:$((b<=PAD_SIZE?1:b-PAD_SIZE))-$((b - 1 + SLICE_SIZE + PAD_SIZE))"

if [[ $REGION_SIZE -eq $SLICE_SIZE ]]
then
  ln --force -s $bamlist $TMPR/bamlist
else
  awk -v FS=/ '{print $NF}' $TMP/local_bamlist | $PARALLEL --jobs=${NUM_THREADS} "$SAMTOOLS index $TMP/bams/{1} && $SAMTOOLS view -b -o $TMPR/bams/{1} $TMP/bams/{1} $PADDED_REGION"
  find $TMPR/bams/ -name "*.bam" | sort > $TMPR/bamlist # Get bamlist for this region
fi

# Increase padded region by read length for graphs
GRAPH_PAD_SIZE=$((PAD_SIZE + 151))
GRAPH_PADDED_REGION="${chrom}:$((b<=GRAPH_PAD_SIZE?1:b-GRAPH_PAD_SIZE))-$((b - 1 + SLICE_SIZE + GRAPH_PAD_SIZE))"

if [[ ! -z  $VCF ]]; then
  $BCFTOOLS view --no-version --output-file $TMPR/region.vcf.gz -Oz $VCF $PADDED_REGION
  $TABIX -f -p vcf $TMPR/region.vcf.gz
fi

echo -n "["

if [[ ! -z $VCF ]] && [[ ${GENOTYPE_ONLY} -ne 0 ]]
then
  ##
  # Iteration 1 (genotype only)
  ##
  mkdir --parents $TMPR/it1

  echo -n "G"

  $GRAPHTYPER construct ${GRAPH}.1 $GENOME --vcf=$VCF ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.1

  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS ${GRAPH}.1 "." \
    --threads=${GRAPHTYPER_THREADS} \
    --no_new_variants \
    --output=$TMPR/it1 \
    --sams=$TMPR/bamlist > ${GT_LOG}.1

  # Create a VCF with all called variants and later join with SV calling
  find $TMPR/it1/ -name "*.hap" > $TMPR/haps1
  $GRAPHTYPER haplotypes ${GRAPH}.1 \
    --haplotypes $TMPR/haps1 \
    --output=$TMP/haps/${chrom}/${region_id}.vcf.gz \
    --skip_breaking_down_extracted_haplotypes \
    --region=$UNPADDED_REGION
  $TABIX -p vcf $TMP/haps/${chrom}/${region_id}.vcf.gz

  hap_calls_vcf=$TMPR/it1/*_calls.vcf.gz

  $GRAPHTYPER vcf_break_down ${GRAPH}.1 ${hap_calls_vcf} \
            --region=$UNPADDED_REGION \
            --output=$TMP/results/${chrom}/${region_id}.vcf.gz

  $TABIX -p vcf $TMP/results/${chrom}/${region_id}.vcf.gz
else
  ##
  # Iteration 1
  ##
  i=1
  mkdir --parents $TMPR/it${i}

  $GRAPHTYPER construct ${GRAPH}.${i} $GENOME ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.${i}
  echo -n "D" # Report discovery iteration

  $GRAPHTYPER discover \
    ${GRAPH}.${i} \
    ${GRAPH_PADDED_REGION} \
    --output=$TMPR/it${i} \
    --sams=$TMPR/bamlist \
    --threads=${GRAPHTYPER_THREADS} \
    --minimum_variant_support=4 \
    --minimum_variant_support_ratio=0.25 > ${GT_LOG}.${i}

  find $TMPR/it${i} -name "*_variant_map" -type f | sort > $TMPR/it${i}/all_variant_maps
  $GRAPHTYPER discovery_vcf ${GRAPH}.${i} $TMPR/it${i}/all_variant_maps --output=$TMPR/discovery.vcf.gz

  if [[ ! -z $VCF ]]; then
    echo -n "V" # Report that the graph was initialized with a VCF file
    $GRAPHTYPER vcf_concatenate --sites_only $TMPR/region.vcf.gz $TMPR/discovery.vcf.gz --output=$TMPR/it${i}/new.vcf.gz
  else
    ln -s $TMPR/discovery.vcf.gz $TMPR/it${i}/new.vcf.gz
  fi

  $TABIX -p vcf $TMPR/it${i}/new.vcf.gz

  ##
  # Iteration 2
  ##
  i=2
  $GRAPHTYPER construct ${GRAPH}.${i} $GENOME --vcf=$TMPR/it1/new.vcf.gz ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.${i}

  echo -n "|D" # Report discovery iteration
  $GRAPHTYPER call \
    $GRAPHTYPER_COMMON_OPTS \
    ${GRAPH}.${i} \
    --threads=${GRAPHTYPER_THREADS} \
    --max_files_open=${MAX_FILES_OPEN} \
    --output=$TMPR/it${i} \
    --sams=$TMPR/bamlist > ${GT_LOG}.${i}

  find $TMPR/it${i}/ -name "*.hap" > $TMPR/haps${i}
  $GRAPHTYPER haplotypes ${GRAPH}.${i} \
    --haplotypes=$TMPR/haps${i} \
    --output=$TMPR/it${i}/haps.vcf.gz

  find $TMPR/it${i}/ -name "*_variant_map" -type f | sort > $TMPR/it${i}/all_variant_maps
  $GRAPHTYPER discovery_vcf ${GRAPH}.${i} $TMPR/it${i}/all_variant_maps --output=$TMPR/it${i}/discovery.vcf.gz
  $GRAPHTYPER vcf_concatenate --sites_only \
    $TMPR/it${i}/haps.vcf.gz \
    $TMPR/it${i}/discovery.vcf.gz \
    --output=$TMPR/it${i}/new.vcf.gz
  $TABIX -p vcf $TMPR/it${i}/new.vcf.gz

  # hash the current variants
  var_hash_curr="$(get_hash $TMPR/it${i}/new.vcf.gz)"

  ##
  # Iteration 3
  # Genotyping-only iteration 1 (cleans graph)
  ##
  i=3
  $GRAPHTYPER construct ${GRAPH}.${i} $GENOME --vcf=$TMPR/it2/new.vcf.gz ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.${i}

  echo -n "|G" # Report genotyping iteration
  $GRAPHTYPER call \
    $GRAPHTYPER_COMMON_OPTS \
    ${GRAPH}.${i} \
    --threads=${GRAPHTYPER_THREADS} \
    --max_files_open=${MAX_FILES_OPEN} \
    --no_new_variants \
    --output=$TMPR/it${i} \
    --sams=$TMPR/bamlist > ${GT_LOG}.${i}

  find $TMPR/it${i}/ -name "*.hap" | sort > $TMPR/haps${i}
  $GRAPHTYPER haplotypes ${GRAPH}.${i} \
    --haplotypes $TMPR/haps${i} \
    --output=$TMPR/it${i}/haps.vcf.gz \
    --skip_breaking_down_extracted_haplotypes

  $TABIX -p vcf $TMPR/it${i}/haps.vcf.gz

  # hash the current variants
  var_hash_prev=${var_hash_curr}
  var_hash_curr="$(get_hash $TMPR/it${i}/haps.vcf.gz)"

  ##
  # Iteration 4
  # Genotyping-only iteration 2 (cleans graph further)
  ##
  if [[ "$var_hash_curr" != "$var_hash_prev" ]]
  then
    i=4
    $GRAPHTYPER construct ${GRAPH}.${i} $GENOME --vcf=$TMPR/it3/haps.vcf.gz ${GRAPH_PADDED_REGION}
    $GRAPHTYPER index ${GRAPH}.${i}
    echo -n "|G" # Report genotyping iteration

    $GRAPHTYPER call \
      ${GRAPHTYPER_COMMON_OPTS} \
      ${GRAPH}.${i} \
      --threads=${GRAPHTYPER_THREADS} \
      --max_files_open=${MAX_FILES_OPEN} \
      --no_new_variants \
      --output=$TMPR/it${i} \
      --sams=$TMPR/bamlist > ${GT_LOG}.${i}

    find $TMPR/it${i}/ -name "*.hap" > $TMPR/haps${i}

    $GRAPHTYPER haplotypes ${GRAPH}.${i} \
      --haplotypes $TMPR/haps${i} \
      --output=$TMPR/it${i}/haps.vcf.gz \
      --skip_breaking_down_extracted_haplotypes

    $TABIX -p vcf $TMPR/it${i}/haps.vcf.gz

    # hash the current variants
    var_hash_prev=${var_hash_curr}
    var_hash_curr="$(get_hash $TMPR/it${i}/haps.vcf.gz)"
  fi

  ##
  # Iteration 5
  # Genotyping-only iteration 3 (make final small variant calls)
  ##
  if [[ "$var_hash_curr" != "$var_hash_prev" ]]
  then
    i=5
    $GRAPHTYPER construct ${GRAPH}.${i} $GENOME --vcf=$TMPR/it4/haps.vcf.gz ${GRAPH_PADDED_REGION}
    $GRAPHTYPER index ${GRAPH}.${i}
    echo -n "|G" # Report genotyping iteration
    $GRAPHTYPER call \
      $GRAPHTYPER_COMMON_OPTS \
      ${GRAPH}.${i} \
      --threads=${GRAPHTYPER_THREADS} \
      --max_files_open=${MAX_FILES_OPEN} \
      --no_new_variants \
      --output=$TMPR/it${i} \
      --sams=$TMPR/bamlist > ${GT_LOG}.${i}

    # Create a VCF with all called variants and later join with SV calling
    find $TMPR/it${i}/ -name "*.hap" > $TMPR/haps${i}

    $GRAPHTYPER haplotypes ${GRAPH}.${i} \
      --haplotypes $TMPR/haps${i} \
      --output=$TMPR/it${i}/haps.vcf.gz \
      --skip_breaking_down_extracted_haplotypes \
      --region=$UNPADDED_REGION

    $TABIX -p vcf $TMPR/it${i}/haps.vcf.gz
  fi

  ##
  # Copy final calls to results
  ##
  cp -p $TMPR/it${i}/haps.vcf.gz $TMP/haps/${chrom}/${region_id}.vcf.gz
  cp -p $TMPR/it${i}/haps.vcf.gz.tbi $TMP/haps/${chrom}/${region_id}.vcf.gz.tbi

  if [[ $(ls $TMPR/it${i}/*_calls.vcf.gz | wc -l) -eq 1 ]]
  then
    ln -s $TMPR/it${i}/*_calls.vcf.gz $TMPR/hap_calls.vcf.gz
  else
    $GRAPHTYPER vcf_merge $TMPR/it${i}/*_calls.vcf.gz --output=$TMPR/hap_calls.vcf.gz
  fi

  $GRAPHTYPER vcf_break_down ${GRAPH}.${i} $TMPR/hap_calls.vcf.gz \
    --region=$UNPADDED_REGION \
    --output=$TMP/results/${chrom}/${region_id}.vcf.gz

  $TABIX -p vcf $TMP/results/${chrom}/${region_id}.vcf.gz
fi

# Clean up
if [[ $CLEAN_UP -ne 0 ]]; then
  rm -r -f $TMPR
fi

echo -n "] "
echo "Completed slice ${chrom}:${b}-$((b - 1 + SLICE_SIZE))"
