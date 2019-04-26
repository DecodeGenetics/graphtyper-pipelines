#!/usr/bin/bash
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
echo -n "["
test `tail -c1 $bamlist` && NUM_SAMPLES=$((NUM_SAMPLES + 1))

# Graphtyper settings
GRAPH=$TMPR/graph
VT_LOG=$TMPR/vt_log
GT_LOG=$TMPR/gt_log

UNPADDED_REGION="${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
PADDED_REGION="${chrom}:$((i<=PAD_SIZE?1:i-PAD_SIZE))-$((i - 1 + SLICE_SIZE + PAD_SIZE))"

if [[ ! -z  $VCF ]]; then
  bcftools view --output-file $TMPR/region.vcf.gz -Oz $VCF $PADDED_REGION
  $TABIX -f $TMPR/region.vcf.gz
fi

while read bamfile; do
  $SAMTOOLS view -b -o $TMPR/bams/$(basename $bamfile) $bamfile $PADDED_REGION
done < $bamlist

find $TMPR/bams/ -name "*.bam" | sort > $TMPR/bamlist # Get bamlist for this region

# Increase padded region by read length for graphs
PAD_SIZE=$((PAD_SIZE + 151))
PADDED_REGION="${chrom}:$((i<=PAD_SIZE?1:i-PAD_SIZE))-$((i - 1 + SLICE_SIZE + PAD_SIZE))"

if [[ ! -z $VCF ]] && [[ ${GENOTYPE_ONLY} -ne 0 ]]
then
  ##
  # Iteration 1 (genotype only)
  ##
  mkdir --parents $TMPR/it1

  echo -n "G"

  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$VCF $PADDED_REGION
  $GRAPHTYPER index $GRAPH

  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
    --threads=${GRAPHTYPER_THREADS}\
    --no_new_variants\
    --output=$TMPR/it1\
    --sams=$TMPR/bamlist > $GT_LOG

  # Create a VCF with all called variants and later join with SV calling
  find $TMPR/it1/ -name "*.hap" > $TMPR/haps1
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes $TMPR/haps1\
    --output=$TMP/haps/${chrom}/${region_id}.vcf.gz\
    --skip_breaking_down_extracted_haplotypes\
    --region=$UNPADDED_REGION
  $TABIX $TMP/haps/${chrom}/${region_id}.vcf.gz

  hap_calls_vcf=$TMPR/it1/*_calls.vcf.gz
  cp ${hap_calls_vcf} $TMP/hap_calls/${chrom}/${region_id}.vcf.gz
  $TABIX $TMP/hap_calls/${chrom}/${region_id}.vcf.gz

  $GRAPHTYPER vcf_break_down $GRAPH ${hap_calls_vcf} \
            --region=$UNPADDED_REGION \
            --output=$TMP/results/${chrom}/${region_id}.vcf.gz

  $TABIX $TMP/results/${chrom}/${region_id}.vcf.gz
else
  ##
  # Iteration 1
  ##
  mkdir --parents $TMPR/it1

  $GRAPHTYPER construct $GRAPH $GENOME $PADDED_REGION
  $GRAPHTYPER index $GRAPH
  echo -n "D" # Report discovery iteration
  $GRAPHTYPER discover $GRAPH $PADDED_REGION \
    --output=$TMPR/it1\
    --sams=$TMPR/bamlist\
    --minimum_variant_support=4\
    --minimum_variant_support_ratio=0.25 > $GT_LOG

  find $TMPR/it1 -name "*_variant_map" -type f | sort > $TMPR/it1/all_variant_maps
  $GRAPHTYPER discovery_vcf $GRAPH $TMPR/it1/all_variant_maps --output=$TMPR/discovery.vcf.gz

  if [[ ! -z $VCF ]]; then
    echo -n "V" # Report that the graph was initialized with a VCF file
    $GRAPHTYPER vcf_concatenate --sites_only $TMPR/region.vcf.gz $TMPR/discovery.vcf.gz --output=$TMPR/it1/new.vcf.gz
  else
    mv $TMPR/discovery.vcf.gz $TMPR/it1/new.vcf.gz
  fi

  $TABIX $TMPR/it1/new.vcf.gz

  # Clear graph
  rm --recursive --force ${GRAPH} ${GRAPH}_gti/

  ##
  # Iteration 2
  ##
  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it1/new.vcf.gz $PADDED_REGION
  $GRAPHTYPER index $GRAPH

  echo -n "|D" # Report discovery iteration
  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
    --threads=${GRAPHTYPER_THREADS}\
    --output=$TMPR/it2\
    --sams=$TMPR/bamlist > $GT_LOG


  find $TMPR/it2/ -name "*.hap" > $TMPR/haps2
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes $TMPR/haps2\
    --output=$TMPR/it2/haps.vcf.gz

  find $TMPR/it2/ -name "*_variant_map" -type f > $TMPR/it2/all_variant_maps
  $GRAPHTYPER discovery_vcf $GRAPH $TMPR/it2/all_variant_maps --output=$TMPR/it2/discovery.vcf.gz
  $GRAPHTYPER vcf_concatenate --sites_only $TMPR/it2/haps.vcf.gz $TMPR/it2/discovery.vcf.gz --output=$TMPR/it2/new.vcf.gz
  $TABIX $TMPR/it2/new.vcf.gz

  # Clear graph and index
  rm -r --force ${GRAPH} ${GRAPH}_gti/

  ##
  # Iteration 3
  # Genotyping-only iteration 1 (cleans graph)
  ##
  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it2/new.vcf.gz $PADDED_REGION
  $GRAPHTYPER index $GRAPH

  echo -n "|G" # Report genotyping iteration
  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
    --threads=${GRAPHTYPER_THREADS}\
    --no_new_variants\
    --output=$TMPR/it3\
    --sams=$TMPR/bamlist > $GT_LOG

  find $TMPR/it3/ -name "*.hap" > $TMPR/haps3
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes $TMPR/haps3\
    --output=$TMPR/it3/haps.vcf.gz\
    --skip_breaking_down_extracted_haplotypes

  $TABIX $TMPR/it3/haps.vcf.gz

  # Clear graph and index
  rm -r --force ${GRAPH} ${GRAPH}_gti/

  ##
  # Iteration 4
  # Genotyping-only iteration 2 (cleans graph further)
  ##
  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it3/haps.vcf.gz $PADDED_REGION
  $GRAPHTYPER index $GRAPH
  echo -n "|G" # Report genotyping iteration
  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
    --threads=${GRAPHTYPER_THREADS}\
    --no_new_variants\
    --output=$TMPR/it4\
    --sams=$TMPR/bamlist > $GT_LOG

  find $TMPR/it4/ -name "*.hap" > $TMPR/haps4
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes $TMPR/haps4\
    --output=$TMPR/it4/haps.vcf.gz\
    --skip_breaking_down_extracted_haplotypes

  $TABIX $TMPR/it4/haps.vcf.gz

  # Clear graph and index
  rm -r --force ${GRAPH} ${GRAPH}_gti/

  ##
  # Iteration 5
  # Genotyping-only iteration 3 (make final small variant calls)
  ##
  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it4/haps.vcf.gz $PADDED_REGION
  $GRAPHTYPER index $GRAPH
  echo -n "|G" # Report genotyping iteration
  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
    --threads=${GRAPHTYPER_THREADS}\
    --no_new_variants\
    --output=$TMPR/it5\
    --sams=$TMPR/bamlist > $GT_LOG

  # Create a VCF with all called variants and later join with SV calling
  find $TMPR/it5/ -name "*.hap" > $TMPR/haps5
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes $TMPR/haps5\
    --output=$TMP/haps/${chrom}/${region_id}.vcf.gz\
    --skip_breaking_down_extracted_haplotypes\
    --region=$UNPADDED_REGION
  $TABIX $TMP/haps/${chrom}/${region_id}.vcf.gz

  hap_calls_vcf=$TMPR/it5/*_calls.vcf.gz
  cp ${hap_calls_vcf} $TMP/hap_calls/${chrom}/${region_id}.vcf.gz
  $TABIX $TMP/hap_calls/${chrom}/${region_id}.vcf.gz

  $GRAPHTYPER vcf_break_down $GRAPH ${hap_calls_vcf} \
            --region=$UNPADDED_REGION \
            --output=$TMP/results/${chrom}/${region_id}.vcf.gz

  $TABIX $TMP/results/${chrom}/${region_id}.vcf.gz
fi

# Clean up
if [[ $CLEAN_UP -ne 0 ]]; then
  rm -r -f $TMPR
fi

echo -n "] "
echo "Completed slice ${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
