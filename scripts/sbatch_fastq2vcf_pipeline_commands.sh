#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/%u/fastq2vcf_pipeline_commands-%j.out
#SBATCH --error=/home/group.kurse/%u/fastq2vcf_pipeline_commands-%j.err

module load python
module load fastqc
module load bwamem2
module load samtools

fastq2vcf_pipeline_commands "$@"
