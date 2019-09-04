#!/usr/bin/env bash

set -e
set -o pipefail

CONFIG=$1
bamlist=$2 # E.g. /nfs/data/.../PN.bam
region=$3 # E.g. chr21:10000

# Parse region
chrom=`echo $region | cut -d":" -f1`
start=`echo $region | cut -d":" -f2`

# Load config
source $CONFIG

end=$((start - 1 + REGION_SIZE))
CHROM_SIZE=`grep -w "^${chrom}" $GENOME.fai | cut -f2`
# Never expand end further than the total length of the chromosome
end=$((end>CHROM_SIZE?CHROM_SIZE:end))

# Pad region
start_padded=$((start - PAD_SIZE < 1 ? 1 : start - PAD_SIZE))
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

Genotyping region:
${UNPADDED_REGION}

Date:
$(date)

Other info:
GENOME=${GENOME}
VCF=${VCF}
BAMSHRINK=${BAMSHRINK}"

mkdir --parents $TMP/results/${chrom} $TMP/haps/${chrom} $TMP/bams

# Clean up temporary directory on failures
trap cleanup 1 2 3 6 15

cleanup()
{
  rm -r --force $TMP
  exit 1
}

# Copy files to local directory to reduce I/O
cp $CONFIG $TMP/config.sh
cp ${TOP_DIR}/call_script.sh $TMP/call_script.sh

if [[ -f "${TOP_DIR}/my_config.sh" ]]
then
  cp ${TOP_DIR}/my_config.sh $TMP/my_config.sh
fi

# Assume it is a single bam file if its extension ends with ".bam"
if [[ "$bamlist" == *.bam ]]
then
  echo "$bamlist" > $TMP/global_bamlist
else
  cp $bamlist $TMP/global_bamlist
fi

# Create local bamlist
cut -f1 $TMP/global_bamlist | awk -v FS=/ -v path="$TMP/bams" '{print path"/"$NF}' > $TMP/local_bamlist

# Copy BAM files to a local directory
if [[ ! -z $BAMSHRINK ]]
then
  echo ${chrom}$'\t'${start}$'\t'${end} > $TMP/region_file

  while IFS=$'\t' read bam cov
  do
    if [[ -z $cov ]]; then
      ${SAMTOOLS} idxstats ${bam} | head -n -1 | awk '{sum+=$3+$4; ref+=$2;} END{print sum/ref}' | xargs echo
    else
      echo ${cov}
    fi
  done < $TMP/global_bamlist > $TMP/global_coverage

  cut -f1 $TMP/global_bamlist > $TMP/global_bamlist2
  mv --force $TMP/global_bamlist2 $TMP/global_bamlist

  $PARALLEL --halt=now,fail=1 --jobs=${NUM_THREADS} --xapply "${BAMSHRINK} {2} {1} 1000 Y 45 {3} {2}.bai $TMP/region_file"\
    :::: $TMP/local_bamlist \
    :::: $TMP/global_bamlist \
    :::: $TMP/global_coverage > $TMP/bamshrink_log_out 2> $TMP/bamshrink_log_err
else
  $PARALLEL -a $TMP/local_bamlist -a <(cut -f1 $TMP/global_bamlist) --halt=now,fail=1 --jobs=${NUM_THREADS} --xapply "${SAMTOOLS} view -o {1} -b {2} ${chrom}:${start_padded}-${end_padded}" > $TMP/samtools_log_out 2> $TMP/samtools_log_err
fi

echo "
Flag meaning:
 - D: Discovery iteration.
 - G: Genotyping-only iteration.
 - V: VCF file was used to initialize the graph.
 - |: Iteration separator."

# Get wall-clock time of copying files to local disk
copy_files_time=$(date +%s)

# Create an array with all regions and paralize the call script
while [[ $start -lt $end ]]
do
  echo "${chrom}:${start}"
  start=$((start + SLICE_SIZE))
done | $PARALLEL --jobs=${NUM_SLICES_RUNNING} --halt=now,fail=1 bash $TMP/call_script.sh $TMP/config.sh $TMP/local_bamlist

# Get wall-clock time of genotyping
genotyping_time=$(date +%s)

# Concatenate all VCF files
$GRAPHTYPER vcf_concatenate $TMP/results/${chrom}/*.vcf.gz --no_sort --output=$TMP/final_small_variants.vcf.gz
$TABIX $TMP/final_small_variants.vcf.gz
$GRAPHTYPER vcf_concatenate $TMP/haps/${chrom}/*.vcf.gz --no_sort --output=$TMP/final_haps.vcf.gz
$TABIX $TMP/final_haps.vcf.gz

# Make sure the output directories exist
mkdir --parents results/${chrom}/ results/.haps/${chrom}

# Move results to output directories
mv $TMP/final_small_variants.vcf.gz results/${chrom}/${region_id}.vcf.gz
mv $TMP/final_small_variants.vcf.gz.tbi results/${chrom}/${region_id}.vcf.gz.tbi
mv $TMP/final_haps.vcf.gz results/.haps/${chrom}/${region_id}.vcf.gz
mv $TMP/final_haps.vcf.gz.tbi results/.haps/${chrom}/${region_id}.vcf.gz.tbi

# Clean up
if [[ $CLEAN_UP -eq 1 ]]; then
  rm -r --force $TMP
fi

echo "
Graphtyper finished genotyping the region successfully. Final results are at:
results/${chrom}/${region_id}.vcf.gz

Total wall-clock time copying files to local disk (seconds):
$((copy_files_time - start_time))

Total wall-clock time of genotyping with Graphtyper (seconds):
$((genotyping_time - copy_files_time))

Total time (seconds):
$((genotyping_time - start_time))
"
