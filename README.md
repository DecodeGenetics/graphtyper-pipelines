# Graphtyper pipelines
This repository has the recommended pipeline script to use for Graphtyper.

The pipeline requires you to have access to the following tools:
 - Graphtyper
 - GNU parallel
 - samtools
 - tabix
 - vt ( https://github.com/atks/vt )

The pipeline script will automatically search for these tools in your PATH directories, but you can also specify their paths by changing the `config.sh` file or create your own `my_config.sh` file. You can also set various parameters in the configuration files.

## Usage
### Short version
```sh
$ bash make_graphtyper_pipeline.sh <BAM> [CONFIG] | bash
```
where BAM can either be a single BAM file or a file with a list of BAM files. CONFIG is the configuration file to use (default: `./config.sh`).

### Long version
#### Running on a computer cluster
The command in the "short version" is a typical use case for someone who wants to run the Graphtyper pipeline a personal computer. If you have a computer cluster, you are likely to be using some workload manager. Integrating the pipeline script with your workload manager should be easy, as each line can be run independently and thus can be run in parallel. For example, if your workload manager of choice is Slurm you can integrate the pipeline with it using:

```sh
bash make_graphtyper_pipeline.sh <BAM> [CONFIG] | \
  while IFS='' read -r line
  do
    srun -c 1 $line &
  done

wait
```

This will run all regions as separate jobs on your cluster using Slurm with one thread allocated.

#### Failed runs
The pipeline automatically detects which VCF files are missing from the results directory. So if you had some failed run you can again run

```sh
$ bash make_graphtyper_pipeline.sh <BAM>
```

to get the commands to run only failed regions. If no commands are outputted, it means all runs have completed.

# License
GNU GPLv3
