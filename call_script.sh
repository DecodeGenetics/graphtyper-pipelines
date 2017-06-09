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
test `tail -c1 NA12878trio` && NUM_SAMPLES=$((NUM_SAMPLES + 1))

# Check which option to use
if [[ $NUM_SAMPLES -le $SMALL_SAMPLE_SIZE ]]
then
  GRAPHTYPER_COMMON_OPTS=$GRAPHTYPER_SMALL_SAMPLE_SIZE_OPTS
else
  GRAPHTYPER_COMMON_OPTS=$GRAPHTYPER_POPULATION_OPTS
fi

# Graphtyper settings
GRAPH=$TMPR/graph
LOGFILE=$TMPR/gt_log
VT_LOG=$TMPR/vt_log

unpadded_region="${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
padded_region="${chrom}:$((i<=PAD_SIZE?1:i-PAD_SIZE))-$((i - 1 + SLICE_SIZE + PAD_SIZE))"

if [[ $INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then
  $VT view -i $padded_region -o $TMPR/region.vcf.gz $VCF
  $TABIX --force $TMPR/region.vcf.gz
fi

for bamfile in `cat $bamlist`
do
  $SAMTOOLS view -b -o $TMPR/bams/$(basename $bamfile) $bamfile $padded_region
done

ls $TMPR/bams/*.bam > $TMPR/bamlist # Get bamlist for this region

# Discovery iteration 1
if [[ INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then
  echo -n "V" # Report that the graph was initialized with a VCF file
  $GRAPHTYPER construct $GRAPH $GENOME $padded_region \
    --vcf=$TMPR/region.vcf.gz\
    --log=$LOGFILE 2>/dev/null
else
  $GRAPHTYPER construct $GRAPH $GENOME --log=$LOGFILE $padded_region 2>/dev/null
fi

$GRAPHTYPER index $GRAPH --log=$LOGFILE
echo -n "D" # Report discovery iteration
$GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
  --output=$TMPR/it1\
  --sams=$TMPR/bamlist\
  --log=$LOGFILE\
  >/dev/null

num_var_before=0

if [[ $INITIALIZE_GRAPH_WITH_VCF -ne 0 ]]; then
  ls $TMPR/it1/*.hap > $TMPR/haps
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes=$TMPR/haps\
    --output=$TMPR/it1/haps.vcf.gz\
    --skip_breaking_down_extracted_haplotypes\
    --log=$LOGFILE

  num_var_before=`zcat $TMPR/region.vcf.gz | grep -v '^#' | wc -l || true`
  $VT cat $TMPR/it1/haps.vcf.gz $TMPR/it1/*_variants.vcf.gz \
    | $VT sort -o $TMPR/new_region_sorted.vcf.gz -
else
  $VT sort -o $TMPR/new_region_sorted.vcf.gz $TMPR/it1/*_variants.vcf.gz
fi

$VT uniq -o $TMPR/new.vcf $TMPR/new_region_sorted.vcf.gz 2> $VT_LOG
num_var_after=`cat $TMPR/new.vcf | grep -v '^#' | wc -l || true`

if [[ $num_var_after -eq 0 ]]
then
  # No variants, just use our current result
  ln -s $TMPR/it1/*_calls.vcf.gz $TMPR/it4/
else
  cat $TMPR/new.vcf | bgzip -c > $TMPR/new_region.vcf.gz
  $TABIX $TMPR/new_region.vcf.gz

  # Clear graph
  rm --force ${GRAPH}
  rm -r --force ${GRAPH}_gti/

  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/new_region.vcf.gz --log=$LOGFILE $padded_region
  $GRAPHTYPER index --log=$LOGFILE $GRAPH

  # Add second discovery iteration when calling a large sample size
  if [[ $NUM_SAMPLES -gt $SMALL_SAMPLE_SIZE ]]; then
    echo -n "|D" # Report discovery iteration
    $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH "." \
      --output=$TMPR/it2\
      --sams=$TMPR/bamlist\
      --log=$LOGFILE\
      >/dev/null

    ls $TMPR/it2/*.hap > $TMPR/haps2
    $GRAPHTYPER haplotypes $GRAPH \
      --haplotypes $TMPR/haps2\
      --output=$TMPR/it2/haps.vcf.gz

    $VT cat $TMPR/it2/haps.vcf.gz $TMPR/it2/*_variants.vcf.gz\
      | $VT sort -o $TMPR/new_region_sorted2.vcf.gz -

    $VT uniq -o $TMPR/new2.vcf $TMPR/new_region_sorted2.vcf.gz 2> $VT_LOG
    cat $TMPR/new2.vcf | bgzip -c > $TMPR/new_region2.vcf.gz
    $TABIX --force $TMPR/new_region2.vcf.gz

    # Clear graph
    rm --force ${GRAPH}
    rm -r --force ${GRAPH}_gti/

    $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/new_region2.vcf.gz --log=$LOGFILE $padded_region
    $GRAPHTYPER index --log=$LOGFILE $GRAPH
  fi

  # Genotyping-only iteration 1 (cleans graph)
  echo -n "|G" # Report genotyping iteration
  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH \
    --no_new_variants\
    --output=$TMPR/it3\
    --sams=$TMPR/bamlist\
    --log=$LOGFILE "." >/dev/null

  ls $TMPR/it3/*.hap > $TMPR/haps3
  $GRAPHTYPER haplotypes $GRAPH \
    --haplotypes $TMPR/haps3\
    --output=$TMPR/it3/haps.vcf.gz\
    --skip_breaking_down_extracted_haplotypes

  $TABIX $TMPR/it3/haps.vcf.gz

  # Clear graph
  rm --force ${GRAPH}
  rm -r --force ${GRAPH}_gti/

  # Genotyping-only iteration 2 (make final calls)
  $GRAPHTYPER construct $GRAPH $GENOME --vcf=$TMPR/it3/haps.vcf.gz --log=$LOGFILE $padded_region
  $GRAPHTYPER index --log=$LOGFILE $GRAPH
  echo -n "|G" # Report genotyping iteration
  $GRAPHTYPER call $GRAPHTYPER_COMMON_OPTS $GRAPH \
    --no_new_variants\
    --output=$TMPR/it4\
    --sams=$TMPR/bamlist\
    --log=$LOGFILE "." >/dev/null
fi

$GRAPHTYPER vcf_break_down $GRAPH $TMPR/it4/*_calls.vcf.gz \
  --region=$unpadded_region \
  --log=$LOGFILE \
  | $VT sort -o $TMP/results/${chrom}/${region_id}.vcf.gz -

$TABIX $TMP/results/${chrom}/${region_id}.vcf.gz

echo -n "] "
echo "Completed slice ${chrom}:${i}-$((i - 1 + SLICE_SIZE))"
