#!/bin/bash

#SBATCH --job-name=TRANSPOSASES                                          # Job name
#SBATCH --ntasks=1                                                          # Run a single task
#SBATCH --cpus-per-task=1                                                   # Number of CPU cores per task
#SBATCH --mem=1G                                                           # Job memory request
#SBATCH --time=3:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=../../logs/116-hpf-transposases/transposases_%j.log       # Standard output and error log
#SBATCH --partition=short

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "" 

# Activate local conda environment
source /home/${USER}/.bashrc

mamba activate snakemake

run_snakemake() {
    local snakefile="$1"
    shift   # remove snakefile from arguments

    if [[ -z "$snakefile" ]]; then
        echo "Usage: run_snakemake <snakefile> [snakemake options]"
        return 1
    fi

    snakemake \
        --snakefile "$snakefile" \
        --profile /scratch/proj/00-default/smk-profiles/short/ \
        --singularity-prefix /scratch/proj/02-compost-microbes/cache/00-singularity \
        --conda-prefix /scratch/proj/02-compost-microbes/cache/00-conda-env \
        --singularity-args "--bind /home:/home,/vast:/vast,/work:/work" \
        --conda-frontend mamba \
        --nolock \
        --rerun-incomplete \
        --keep-incomplete \
        --allow-ambiguity \
        --keep-going \
        --rerun-triggers mtime \
        --jobs 600 \
        "$@"
}

run_snakemake main.smk --forcerun blast_makedb