#!/usr/bin/bash
set -e
set -o pipefail

CONFIG=$1
bam=$2 # E.g. /nfs/data/.../PN.bam
bambai=${bam}.bai
region=$3 # E.g. chr21:10000

# Parse region
chrom=`echo $region | cut -d":" -f1`
start=`echo $region | cut -d":" -f2`
i=$start

# Load config
source $CONFIG

end=$((start - 1 + REGION_SIZE))
CHROM_SIZE=`grep -w "^${chrom}" $GENOME.fai | cut -f2`
# Never expand end further than the total length of the chromosome
end=$((end>CHROM_SIZE?CHROM_SIZE:end))

# Pad region
start_padded=$((i - PAD_SIZE < 1 ? 1 : i - PAD_SIZE))
end_padded=$((end + PAD_SIZE > CHROM_SIZE ? CHROM_SIZE : end + PAD_SIZE))

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
${chrom}:${start}-${end}

Date:
$(date)

Start time (UNIX seconds):
$start_time"

mkdir -p $TMP/results/${chrom}
mkdir -p $TMP/bams

# Clean up temporary directory
if [[ $CLEAN_UP -ne 0 ]]; then
  trap "rm -r --force $TMP; exit " 0 1 2 15
fi

# Copy files to local directory to reduce I/O
cp $CONFIG $TMP/config.sh
cp ${TOP_DIR}/call_script.sh $TMP/call_script.sh

if [[ -f "${TOP_DIR}/my_config.sh" ]]
then
  cp ${TOP_DIR}/my_config.sh $TMP/my_config.sh
fi

# Assume it is a single bam file if its extension ends with ".bam"
if [[ "$bam" == *.bam ]]
then
  echo "$bam" > $TMP/global_bamlist
else
  cp $bam $TMP/global_bamlist
fi

for bamfile in `cat $TMP/global_bamlist`
do
  echo $TMP/bams/$(basename $bamfile)
done > $TMP/local_bamlist

# Copy BAM files to a local directory
$PARALLEL --jobs=${NUM_THREADS} --xapply $SAMTOOLS\
  ::: view\
  ::: -o\
  :::: $TMP/local_bamlist\
  ::: -b\
  :::: $TMP/global_bamlist\
  ::: "${chrom}:${start_padded}-${end_padded}"


# Index BAM files
$PARALLEL --jobs=${NUM_THREADS} $SAMTOOLS\
  ::: index\
  :::: $TMP/local_bamlist

# Get wall-clock time of copying files to local disk
copy_files_time=$(date +%s)

echo "
Total wall-clock time copying files to local disk (seconds):
$((copy_files_time - start_time))"

# Create an array with all regions
regions=()

while [[ i -lt $end ]]
do
  regions+=("${chrom}:${i}")
  i=$((i + SLICE_SIZE))
done

echo "
Flag meaning:
 - D: Discovery iteration.
 - G: Genotyping-only iteration.
 - V: VCF file was used to initialize the graph.
 - |: Iteration separator."

# Paralize the call script
$PARALLEL --jobs=$NUM_SLICES_RUNNING --halt=now,fail=1 bash\
  ::: $TMP/call_script.sh\
  ::: $TMP/config.sh\
  ::: $TMP/local_bamlist\
  ::: $(echo ${regions[*]})

# Get wall-clock time of genotyping
genotyping_time=$(date +%s)

echo "
Total wall-clock time of genotyping with Graphtyper (seconds):
$((genotyping_time - copy_files_time))

Total time (seconds):
$((genotyping_time - start_time))
"

# Concatenate all VCF files
find $TMP/results/${chrom}/ -name "*.vcf.gz" -type f | sort > $TMP/output_files

$VT cat -L $TMP/output_files \
  | $VT sort -o $TMP/final.vcf.gz -m local -w 200 -

$TABIX $TMP/final.vcf.gz

# Move final results
mkdir -p ${RESULTS}/${chrom}
mv $TMP/final.vcf.gz ${RESULTS}/${chrom}/${region_id}.vcf.gz.tmp
mv $TMP/final.vcf.gz.tbi ${RESULTS}/${chrom}/${region_id}.vcf.gz.tbi.tmp

# Check if this region was already been genotyped
if [[ -f "${RESULTS}/${chrom}/${region_id}.vcf.gz" ]]
then
  echo "
WARNING: The output file already exists '${RESULTS}/${chrom}/${region_id}.vcf.gz'.
         The old file will be moved to '${RESULTS}/${chrom}/${region_id}.vcf.gz.bak'.
"
  mv --force ${RESULTS}/${chrom}/${region_id}.vcf.gz ${RESULTS}/${chrom}/${region_id}.vcf.gz.bak
  mv --force ${RESULTS}/${chrom}/${region_id}.vcf.gz.tbi ${RESULTS}/${chrom}/${region_id}.vcf.gz.tbi.bak
fi


mv ${RESULTS}/${chrom}/${region_id}.vcf.gz.tmp ${RESULTS}/${chrom}/${region_id}.vcf.gz
mv ${RESULTS}/${chrom}/${region_id}.vcf.gz.tbi.tmp ${RESULTS}/${chrom}/${region_id}.vcf.gz.tbi

echo "
Graphtyper finished genotyping the region successfully. Final results are at:
${RESULTS}/${chrom}/${region_id}.vcf.gz
"
