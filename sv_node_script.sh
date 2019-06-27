#!/usr/bin/env bash

set -e
set -o pipefail

CONFIG=$1
bam=$2 # E.g. /nfs/data/.../PN.bam or a bamlist
region=$3 # E.g. chr21:10000

# Parse region
chrom=`echo $region | cut -d":" -f1`
start=`echo $region | cut -d":" -f2`

# Load config
source $CONFIG

if [[ ! -f $SV_VCF ]]; then
  echo "ERROR: Failed to find SV_VCF=${SV_VCF}" >&2
  exit 1
fi

SV_REGION_SIZE=1000000 # 1MB
end=$((start - 1 + SV_REGION_SIZE))
CHROM_SIZE=`grep -w "^${chrom}" $GENOME.fai | cut -f2`
# Never expand end further than the total length of the chromosome
end=$((end>CHROM_SIZE?CHROM_SIZE:end))

if [[ ${SV_RESULTS_DIR_SUFFIX_SAMPLE1} -ne 0 ]]; then
  SUFFIX=`head -n 1 ${bam} | cut -f1 | xargs samtools view -H | grep -oP "SM:[A-Z0-9]+" | head -n 1 | cut -d":" -f2`
  SV_RESULTS_DIR="${SV_RESULTS_DIR}${SUFFIX}"
fi

# Pad region
start_padded=$((start - PAD_SIZE < 1 ? 1 : start - PAD_SIZE))
PAD_SIZE=$((PAD_SIZE + 200000)) # Extra padding at end in SV genotyping
end_padded=$((end + PAD_SIZE > CHROM_SIZE ? CHROM_SIZE : end + PAD_SIZE))

UNPADDED_REGION="${chrom}:${start}-${end}"
PADDED_REGION="${chrom}:${start_padded}-${end_padded}"

# Get ID for the region
region_id=`printf "%09d" $start`"-"`printf "%09d" $end`

# Get timing of start
start_time=$(date +%s)

# Create a temporary directories
TMP=$(mktemp --directory ${TMP_FORMAT})
echo "Hostname:
$(hostname | cut -d"." -f1)

Temporary directory:
$TMP

Graphtyper binary:
$(realpath ${GRAPHTYPER})

SV genotyping region:
${UNPADDED_REGION}

SV padded region:
${PADDED_REGION}

SV results dir:
${SV_RESULTS_DIR}

Date:
$(date)
"

mkdir --parents $TMP/results/${chrom} $TMP/bams

trap cleanup 1 2 3 6 15

cleanup()
{
  rm -r --force $TMP
  exit 1
}

# Assume it is a single bam file if its extension ends with ".bam"
if [[ "$bam" == *.bam ]]
then
  echo "$bam" > $TMP/global_bamlist
else
  cp $bam $TMP/global_bamlist
fi

while IFS=$'\t' read bamfile cov; do
  echo $TMP/bams/$(basename $bamfile)
done < $TMP/global_bamlist > $TMP/local_bamlist

bamlist_size=`cat $TMP/local_bamlist | wc -l`

echo `date`" INFO: Bamlist size is ${bamlist_size}"

# Make sure the results directory exists
mkdir --parents ${SV_RESULTS_DIR}/${chrom}/

# SV graph construction
if [[ -f haps/${chrom}.vcf.gz ]]
then
  cat <($VT view -H $SV_VCF | head -n -1) <($VT view -H haps/${chrom}.vcf.gz | grep -P "^##INFO") <($VT view -H $SV_VCF | tail -n 1) > $TMP/vcf_header
  zcat $SV_VCF haps/${chrom}.vcf.gz | awk -v chrom=${chrom} -v start=$start -v end=$end 'substr($1,1,1) != "#" && $1 == chrom && $2 >= start && $2 <= end' | sort -k2,2n | cat $TMP/vcf_header - > $TMP/input_with_svs.vcf
else
  zcat $SV_VCF | grep "^#" > $TMP/vcf_header
  zcat $SV_VCF | awk -v chrom=${chrom} -v start=$start -v end=$end 'substr($1,1,1) != "#" && $1 == chrom && $2 >= start && $2 <= end' | sort -k2,2n | cat $TMP/vcf_header - > $TMP/input_with_svs.vcf
fi

cat $TMP/input_with_svs.vcf | bgzip -c > $TMP/input_with_svs.vcf.gz
$TABIX $TMP/input_with_svs.vcf.gz
GRAPH="$TMP/input_with_svs"

echo `date`" INFO: Constructing graph ${GRAPH}"
$GRAPHTYPER construct $GRAPH $GENOME $PADDED_REGION --vcf=$TMP/input_with_svs.vcf.gz --sv_graph
echo `date`" INFO: Indexing graph ${GRAPH}"
$GRAPHTYPER index $GRAPH
echo `date`" INFO: Calling graph ${GRAPH}"

# Copy BAM files to a local directory
if [[ ! -z $BAMSHRINK ]]; then
  echo ${chrom}$'\t'${start}$'\t'$((end + 50000)) > $TMP/region_file

  while IFS=$'\t' read bam cov; do
    if [[ -z $cov ]]; then
      echo "WARNING: Found a BAM with no avgCovByReadLen calculated. Pre-calculating it is more efficient." >&2
      cov=`$SAMTOOLS idxstats ${bam} | head -n -1 | awk '{sum+=$3+$4; ref+=$2;} END{print sum/ref}'`
    fi

    echo ${cov}
  done < $TMP/global_bamlist > $TMP/global_coverage

  cut -f1 $TMP/global_bamlist > $TMP/global_bamlist2
  mv --force $TMP/global_bamlist2 $TMP/global_bamlist

  $PARALLEL --halt=now,fail=1 --jobs=${NUM_THREADS_SV_CALLING} --xapply "${BAMSHRINK} {2} {1} 1000 Y 45 {3} {2}.bai $TMP/region_file && $GRAPHTYPER call $GRAPH \".\" --threads=1 --no_new_variants --one_genotype_per_haplotype --output=$TMP/sv_results --sam={1} && rm {1}"\
    :::: $TMP/local_bamlist\
    :::: $TMP/global_bamlist\
    :::: $TMP/global_coverage
else
  # TOmaybeDO: support samtools view as well
  echo "ERROR: bamShrink required for SV calling" >&2
  exit 1
fi

# Merge VCF files
find $TMP/sv_results/ -maxdepth 1 -name "*_calls.vcf.gz" -type f | sort > $TMP/all_calls
#split --lines=1000 $TMP/all_calls $TMP/part_
#find $TMP/ -maxdepth 1 -name "part_*" -type f | sort > $TMP/all_parts
#$PARALLEL --halt=now,fail=1 --jobs=1 -a $TMP/all_parts "$GRAPHTYPER vcf_merge --file_list={1} --output={1}.vcf.gz"
#find $TMP/ -maxdepth 1 -name "part_*.vcf.gz" -type f | sort > $TMP/all_vcfs
$GRAPHTYPER vcf_merge --file_list=$TMP/all_calls --output=$TMP/graphtyper_calls_merged.vcf.gz

# Remove variants outside of the genotyping region
zcat $TMP/graphtyper_calls_merged.vcf.gz | awk -v start=$start -v end=$end 'substr($1,1,1) == "#" || ($2 >= start && $2 <= end)' | bgzip -c > $TMP/graphtyper_calls_merged_in_region.vcf.gz
$TABIX -p vcf -f $TMP/graphtyper_calls_merged_in_region.vcf.gz

mv $TMP/graphtyper_calls_merged_in_region.vcf.gz ${SV_RESULTS_DIR}/${chrom}/${region_id}.raw.vcf.gz
mv $TMP/graphtyper_calls_merged_in_region.vcf.gz.tbi ${SV_RESULTS_DIR}/${chrom}/${region_id}.raw.vcf.gz.tbi

# Clean up
if [[ $CLEAN_UP -eq 1 ]]; then
  rm -r --force $TMP
fi

end_time=$(date +%s)

echo "
Graphtyper finished SV genotyping the region successfully. Final results are at:
${SV_RESULTS_DIR}/${chrom}/${region_id}.vcf.gz

Total time (seconds):
$((end_time - start_time))
"
