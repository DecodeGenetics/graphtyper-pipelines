## Test dataset of GraphTyper SV genotyping

```sh
# Construct a graph containing some SVs defined in svs.vcf.gz on a contig called "chr20"
graphtyper construct my_graph reference.fa --vcf=svs.vcf.gz --sv_graph chr20

# Build an index
graphtyper index my_graph

# Genotype SVs on 4 simulated samples (50x coverage) and put them in a directory named "calls"
mkdir -p calls
graphtyper call my_graph . --sam=SAMP1.bam --threads=4 --no_new_variants --output=calls
graphtyper call my_graph . --sam=SAMP2.bam --threads=4 --no_new_variants --output=calls
graphtyper call my_graph . --sam=SAMP3.bam --threads=4 --no_new_variants --output=calls
graphtyper call my_graph . --sam=SAMP4.bam --threads=4 --no_new_variants --output=calls

# Merge results into a single VCF and update INFO fields
graphtyper vcf_merge calls/*.vcf.gz --output=graphtyper.final.vcf.gz

# Extract the aggregated SV genotyping model and variants that passed filters
zgrep "^#\|SVMODEL=AGGREGATED" graphtyper.final.vcf.gz | grep -P "^#|\tPASS\t" > graphtyper.final.agg.PASS.vcf

# Compare to expected output (assuming using graphtyper version with SHA1 hash b9c603e). No output means the files are the same.
diff graphtyper.final.agg.PASS.vcf expected_b9c603e.vcf
```
