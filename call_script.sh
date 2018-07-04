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

mkdir -p $TMPR/bams
mkdir -p $TMPR/it1
mkdir -p $TMPR/it2
mkdir -p $TMPR/it3
mkdir -p $TMPR/it4

# Clean up after calling
if [[ $CLEAN_UP -ne 0 ]]; then
  trap "rm -r -f $TMPR; exit " 0 1 2 15
fi

NUM_SAMPLES=`cat $bamlist | wc -l`
echo -n "["
test `tail -c1 $bamlist` && NUM_SAMPLES=$((NUM_SAMPLES + 1))

# Check which option to use
if [[ $NUM_SAMPLES -le $SMALL_SAMPLE_SIZE ]]
then
  GRAPHTYPER_COMMON_OPTS=$GRAPHTYPER_SMALL_SAMPLE_SIZE_OPTS
else
  GRAPHTYPER_COMMON_OPTS=$GRAPHTYPER_POPULATION_OPTS
fi

# Graphtyper settings
GRAPH=$TMPR/graph
VT_LOG=$TMPR/vt_log
GT_LOG=$TMPR/gt_log

unpadded_region="${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
padded_region="${chrom}:$((i<=PAD_SIZE?1:i-PAD_SIZE))-$((i - 1 + SLICE_SIZE + PAD_SIZE))"

if [[ $INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then
  bcftools view --output-file $TMPR/region.vcf.gz -Oz $VCF $padded_region
  $TABIX -f $TMPR/region.vcf.gz
fi

for bamfile in `cat $bamlist`
do
  $SAMTOOLS view -b -o $TMPR/bams/$(basename $bamfile) $bamfile $padded_region
done

find $TMPR/bams/ -name "*.bam" | sort > $TMPR/bamlist # Get bamlist for this region

##
# Iteration 1
##
if [[ $INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then
  echo -n "V" # Report that the graph was initialized with a VCF file
  $GRAPHTYPER construct $GRAPH $GENOME $padded_region --vcf=$TMPR/region.vcf.gz
else
  $GRAPHTYPER construct $GRAPH $GENOME $padded_region
fi

mkdir --parents $TMPR/it1

$GRAPHTYPER index $GRAPH
echo -n "D" # Report discovery iteration
$GRAPHTYPER discover $GRAPH $padded_region \
  --output=$TMPR/it1\
  --sams=$TMPR/bamlist > $GT_LOG

find $TMPR/it1 -name "*_variant_map" -type f | sort > $TMPR/it1/all_variant_maps
$GRAPHTYPER discovery_vcf $GRAPH $TMPR/it1/all_variant_maps --output=$TMPR/discovery.vcf.gz

num_var_before=0

if [[ $INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then
  num_var_before=`zcat $TMPR/region.vcf.gz | grep -v '^#' | wc -l || true`
  $GRAPHTYPER vcf_concatenate $TMPR/region.vcf.gz $TMPR/discovery.vcf.gz --output=$TMPR/it1/new.vcf.gz
else
  mv $TMPR/discovery.vcf.gz $TMPR/it1/new.vcf.gz
fi

$TABIX $TMPR/it1/new.vcf.gz

num_var_after=`zgrep -v "^#" $TMPR/it1/new.vcf.gz | wc -l || true`

# Clear graph
rm --force ${GRAPH}
rm -r --force ${GRAPH}_gti/

##
# Iteration 2
##
$GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it1/new.vcf.gz $padded_region
$GRAPHTYPER index $GRAPH

# Add second discovery iteration when calling a large sample size
echo -n "|D" # Report discovery iteration
$GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
  --output=$TMPR/it2\
  --sams=$TMPR/bamlist > $GT_LOG

find $TMPR/it2/ -name "*.hap" | sort > $TMPR/haps2
$GRAPHTYPER haplotypes $GRAPH \
  --haplotypes $TMPR/haps2\
  --output=$TMPR/it2/haps.vcf.gz

find $TMPR/it2/ -name "*_variant_map" -type f | sort > $TMPR/it2/all_variant_maps
$GRAPHTYPER discovery_vcf $GRAPH $TMPR/it2/all_variant_maps --output=$TMPR/it2/discovery.vcf.gz
$GRAPHTYPER vcf_concatenate $TMPR/it2/haps.vcf.gz $TMPR/it2/discovery.vcf.gz --output=$TMPR/it2/new.vcf.gz
$TABIX $TMPR/it2/new.vcf.gz

# Clear graph
rm --recursive --force ${GRAPH}
rm -r --force ${GRAPH}_gti/

##
# Iteration 3
# Genotyping-only iteration 1 (cleans graph)
##
$GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it2/new.vcf.gz $padded_region
$GRAPHTYPER index $GRAPH

echo -n "|G" # Report genotyping iteration
$GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
  --no_new_variants\
  --output=$TMPR/it3\
  --sams=$TMPR/bamlist > $GT_LOG

find $TMPR/it3/ -name "*.hap" | sort > $TMPR/haps3
$GRAPHTYPER haplotypes $GRAPH \
  --haplotypes $TMPR/haps3\
  --output=$TMPR/it3/haps.vcf.gz\
  --skip_breaking_down_extracted_haplotypes

$TABIX $TMPR/it3/haps.vcf.gz

# Clear graph
rm --force ${GRAPH}
rm -r --force ${GRAPH}_gti/

# Genotyping-only iteration 2 (make final calls)
$GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it3/haps.vcf.gz $padded_region
$GRAPHTYPER index $GRAPH
echo -n "|G" # Report genotyping iteration
$GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH \
  --no_new_variants\
  --output=$TMPR/it4\
  --sams=$TMPR/bamlist\
  "."  > $GT_LOG


$GRAPHTYPER vcf_break_down $GRAPH $TMPR/it4/*_calls.vcf.gz \
            --region=$unpadded_region \
            --output=$TMP/results/${chrom}/${region_id}.vcf.gz

$TABIX $TMP/results/${chrom}/${region_id}.vcf.gz

echo -n "] "
echo "Completed slice ${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
