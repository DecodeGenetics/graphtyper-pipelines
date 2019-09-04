#!/usr/bin/env bash

set -e
set -o pipefail

CONFIG=$1
bamlist=$2
region=$3 # E.g. chr21:1
TMP=`dirname $bamlist`

source $CONFIG

chrom=`echo $region | cut -d":" -f1`
i=`echo $region | cut -d":" -f2`
region_id=`printf "%09d" $i`"-"`printf "%09d" $((i - 1 + SLICE_SIZE))`

TMPR="$TMP/${chrom}/${region_id}"

mkdir --parents $TMPR/bams $TMPR/it{1,2,3,4,5}

# Clean up after calling
if [[ $CLEAN_UP -ne 0 ]]; then
  trap "rm -r -f $TMPR; exit 1" 1 2 15
fi

NUM_SAMPLES=`cat $bamlist | wc -l`
test `tail -c1 $bamlist` && NUM_SAMPLES=$((NUM_SAMPLES + 1))
NUM_POOLS="$((NUM_POOLS<=NUM_SAMPLES?NUM_POOLS:NUM_SAMPLES))"

# Graphtyper settings
GRAPH=$TMPR/graph
GT_LOG=$TMPR/gt_log

UNPADDED_REGION="${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
PADDED_REGION="${chrom}:$((i<=PAD_SIZE?1:i-PAD_SIZE))-$((i - 1 + SLICE_SIZE + PAD_SIZE))"

if [[ $REGION_SIZE -eq $SLICE_SIZE ]]
then
  ln -s $bamlist $TMPR/bamlist
else
  awk -v FS=/ '{print $NF}' $TMP/local_bamlist | $PARALLEL --jobs=${NUM_THREADS} "$SAMTOOLS index $TMP/bams/{1} && $SAMTOOLS view -b -o $TMPR/bams/{1} $TMP/bams/{1} $PADDED_REGION"
  find $TMPR/bams/ -name "*.bam" | sort > $TMPR/bamlist # Get bamlist for this region
fi

# Split bamlist into pools if NUM_POOLS is more than one
if [[ $NUM_POOLS -le 1 ]]
then
  echo "$TMPR/bamlist" > $TMPR/pools
else
  mkdir --parents $TMPR/bamlists/
  split --number="l/${NUM_POOLS}" --suffix-length=3 $TMPR/bamlist $TMPR/bamlists/p.
  find $TMPR/bamlists -type f | sort > $TMPR/pools
fi

# Increase padded region by read length for graphs
GRAPH_PAD_SIZE=$((PAD_SIZE + 151))
GRAPH_PADDED_REGION="${chrom}:$((i<=GRAPH_PAD_SIZE?1:i-GRAPH_PAD_SIZE))-$((i - 1 + SLICE_SIZE + GRAPH_PAD_SIZE))"

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
    --threads=${GRAPHTYPER_THREADS}\
    --no_new_variants\
    --output=$TMPR/it1\
    --sams=$TMPR/bamlist > ${GT_LOG}.1

  # Create a VCF with all called variants and later join with SV calling
  find $TMPR/it1/ -name "*.hap" > $TMPR/haps1
  $GRAPHTYPER haplotypes ${GRAPH}.1 \
    --haplotypes $TMPR/haps1\
    --output=$TMP/haps/${chrom}/${region_id}.vcf.gz\
    --skip_breaking_down_extracted_haplotypes\
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
  mkdir --parents $TMPR/it1

  $GRAPHTYPER construct ${GRAPH}.1 $GENOME ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.1
  echo -n "D" # Report discovery iteration

  $PARALLEL --halt=now,fail=1 --jobs=${NUM_POOL_THREADS} --arg-file=$TMPR/pools "$GRAPHTYPER discover \
    ${GRAPH}.1 \
    ${GRAPH_PADDED_REGION} \
    --output=$TMPR/it1 \
    --sams={1} \
    --minimum_variant_support=4 \
    --minimum_variant_support_ratio=0.25" > ${GT_LOG}.1

  find $TMPR/it1 -name "*_variant_map" -type f | sort > $TMPR/it1/all_variant_maps
  $GRAPHTYPER discovery_vcf ${GRAPH}.1 $TMPR/it1/all_variant_maps --output=$TMPR/discovery.vcf.gz

  if [[ ! -z $VCF ]]; then
    echo -n "V" # Report that the graph was initialized with a VCF file
    $GRAPHTYPER vcf_concatenate --sites_only $TMPR/region.vcf.gz $TMPR/discovery.vcf.gz --output=$TMPR/it1/new.vcf.gz
  else
    ln -s $TMPR/discovery.vcf.gz $TMPR/it1/new.vcf.gz
  fi

  $TABIX -p vcf $TMPR/it1/new.vcf.gz

  ##
  # Iteration 2
  ##
  $GRAPHTYPER construct ${GRAPH}.2 $GENOME --vcf=$TMPR/it1/new.vcf.gz ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.2

  echo -n "|D" # Report discovery iteration
  $PARALLEL --halt=now,fail=1 --jobs=${NUM_POOL_THREADS} --arg-file=$TMPR/pools "$GRAPHTYPER call \
    $GRAPHTYPER_COMMON_OPTS \
    ${GRAPH}.2 \
    . \
    --threads=${GRAPHTYPER_THREADS} \
    --output=$TMPR/it2 \
    --sams={1}" > ${GT_LOG}.2

  find $TMPR/it2/ -name "*.hap" > $TMPR/haps2
  $GRAPHTYPER haplotypes ${GRAPH}.2 \
    --haplotypes=$TMPR/haps2\
    --output=$TMPR/it2/haps.vcf.gz

  find $TMPR/it2/ -name "*_variant_map" -type f > $TMPR/it2/all_variant_maps
  $GRAPHTYPER discovery_vcf ${GRAPH}.2 $TMPR/it2/all_variant_maps --output=$TMPR/it2/discovery.vcf.gz
  $GRAPHTYPER vcf_concatenate --sites_only $TMPR/it2/haps.vcf.gz $TMPR/it2/discovery.vcf.gz --output=$TMPR/it2/new.vcf.gz
  $TABIX -p vcf $TMPR/it2/new.vcf.gz

  ##
  # Iteration 3
  # Genotyping-only iteration 1 (cleans graph)
  ##
  $GRAPHTYPER construct ${GRAPH}.3 $GENOME --vcf=$TMPR/it2/new.vcf.gz ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.3

  echo -n "|G" # Report genotyping iteration
  $PARALLEL --halt=now,fail=1 --jobs=${NUM_POOL_THREADS} --arg-file=$TMPR/pools "$GRAPHTYPER call \
    $GRAPHTYPER_COMMON_OPTS \
    ${GRAPH}.3 \
    . \
    --threads=${GRAPHTYPER_THREADS} \
    --no_new_variants \
    --output=$TMPR/it3 \
    --sams={1}" > ${GT_LOG}.3

  find $TMPR/it3/ -name "*.hap" > $TMPR/haps3
  $GRAPHTYPER haplotypes ${GRAPH}.3 \
    --haplotypes $TMPR/haps3\
    --output=$TMPR/it3/haps.vcf.gz\
    --skip_breaking_down_extracted_haplotypes

  $TABIX -p vcf $TMPR/it3/haps.vcf.gz

  ##
  # Iteration 4
  # Genotyping-only iteration 2 (cleans graph further)
  ##
  $GRAPHTYPER construct ${GRAPH}.4 $GENOME --vcf=$TMPR/it3/haps.vcf.gz ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.4
  echo -n "|G" # Report genotyping iteration

  $PARALLEL --halt=now,fail=1 --jobs=${NUM_POOL_THREADS} --arg-file=$TMPR/pools "$GRAPHTYPER call \
    ${GRAPHTYPER_COMMON_OPTS} \
    ${GRAPH}.4 \
    . \
    --threads=${GRAPHTYPER_THREADS} \
    --no_new_variants \
    --output=$TMPR/it4 \
    --sams={1}" > ${GT_LOG}.4

  find $TMPR/it4/ -name "*.hap" > $TMPR/haps4

  $GRAPHTYPER haplotypes ${GRAPH}.4 \
    --haplotypes $TMPR/haps4\
    --output=$TMPR/it4/haps.vcf.gz\
    --skip_breaking_down_extracted_haplotypes

  $TABIX -p vcf $TMPR/it4/haps.vcf.gz

  ##
  # Iteration 5
  # Genotyping-only iteration 3 (make final small variant calls)
  ##
  $GRAPHTYPER construct ${GRAPH}.5 $GENOME --vcf=$TMPR/it4/haps.vcf.gz ${GRAPH_PADDED_REGION}
  $GRAPHTYPER index ${GRAPH}.5
  echo -n "|G" # Report genotyping iteration
  $PARALLEL --halt=now,fail=1 --jobs=${NUM_POOL_THREADS} --arg-file=$TMPR/pools "$GRAPHTYPER call \
    $GRAPHTYPER_COMMON_OPTS \
    ${GRAPH}.5 \
    . \
    --threads=${GRAPHTYPER_THREADS} \
    --no_new_variants \
    --output=$TMPR/it5 \
    --sams={1}" > ${GT_LOG}.5

  # Create a VCF with all called variants and later join with SV calling
  find $TMPR/it5/ -name "*.hap" > $TMPR/haps5

  $GRAPHTYPER haplotypes ${GRAPH}.5 \
    --haplotypes $TMPR/haps5\
    --output=$TMP/haps/${chrom}/${region_id}.vcf.gz\
    --skip_breaking_down_extracted_haplotypes\
    --region=$UNPADDED_REGION

  $TABIX -p vcf $TMP/haps/${chrom}/${region_id}.vcf.gz

  if [[ $NUM_POOLS -le 1 ]]
  then
    ln -s $TMPR/it5/*_calls.vcf.gz $TMPR/hap_calls.vcf.gz
  else
    $GRAPHTYPER vcf_merge $TMPR/it5/*_calls.vcf.gz --output=$TMPR/hap_calls.vcf.gz
  fi

  $GRAPHTYPER vcf_break_down ${GRAPH}.5 $TMPR/hap_calls.vcf.gz \
    --region=$UNPADDED_REGION \
    --output=$TMP/results/${chrom}/${region_id}.vcf.gz

  $TABIX -p vcf $TMP/results/${chrom}/${region_id}.vcf.gz
fi

# Clean up
if [[ $CLEAN_UP -ne 0 ]]; then
  rm -r -f $TMPR
fi

echo -n "] "
echo "Completed slice ${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
