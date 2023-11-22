#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/%u/fastq2vcf_pipeline_commands-%j.out
#SBATCH --error=/home/group.kurse/%u/fastq2vcf_pipeline_commands-%j.err

cd $HOME

module load python/3.8.7
module load fastqc
module load bwamem2/2.2.1
module load samtools/1.13

fastq2vcf_pipeline_commands $@
