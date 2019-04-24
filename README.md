# Graphtyper pipelines
This repository has the recommended pipeline scripts to use for Graphtyper. The scripts depends on the following tools:
 - Graphtyper
 - GNU parallel
 - samtools
 - tabix

The pipeline scripts will automatically search for these tools in your PATH directories, but you can also specify their paths by changing the `config.sh` file or create your own `my_config.sh` file. You can also set various variables/parameters in the configuration files (which are documented there). The only required variable is the `GENOME` variable, which you must set to your reference genome FASTA file.

## Usage
### Short version
Small variants genotyping
```sh
bash make_graphtyper_pipeline.sh <BAM> [CONFIG] | bash
```
where BAM can either be a single BAM file or a file with a list of BAM files. CONFIG is the configuration file to use (default: `./config.sh`). The SV genotyping is run similarly:
```sh
bash make_graphtyper_sv_pipeline.sh <BAM> [CONFIG]
```

### Long version
#### Running on a computer cluster
The command in the "short version" is a typical use case for someone who wants to run the Graphtyper pipelines on a personal computer. If you have a computer cluster, you are likely to be using some workload manager. Integrating the pipeline script with your workload manager should be easy, as each line can be run independently and thus can be run in parallel. For example, if your workload manager of choice is Slurm you could use:
```sh
bash make_graphtyper_[sv_]pipeline.sh <BAM> [CONFIG] | \
  while IFS='' read -r line
  do
    srun -c 4 $line &
  done

wait
```

This will run all regions as separate jobs on your cluster using Slurm with four threads allocated.

#### Failed runs
The pipeline automatically detects which VCF files are missing from the results directory. So if you had some failed runs you can again run:
```sh
bash make_graphtyper_[sv_]pipeline.sh <BAM>
```

to get the commands to run only failed regions. If no commands are outputted, it means all runs have completed. If your runs fail becaue of memory, we recommend decreasing the `NUM_SLICES_RUNNING` variable to reduce the memory requirement without having any effect on quality.

## Example
A test example config can be found in `test_config.sh`. You can run the test with:

```sh
make -f pipelines.make test
```

If everything is configured properly the test results should appear in `test/results/20/000000500-000009499.vcf.gz`. That output should be the same as in `test/expected.vcf.gz`.

# License
GNU GPLv3
