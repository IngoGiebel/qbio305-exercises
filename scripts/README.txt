################################################################################
### Required scripts and parameter files
### Must be executable and located within the $PATH environment.
################################################################################

Shell scripts (within the scripts subdirectory):
- sbatch_fastq2vcf_pipeline_commands.sh
- fastq2vcf_pipeline_commands

Parameter script (within the scripts subdirectory):
- fastq2vcf_pipeline_commands.ini

Input groups file (required for the qualimap multi-bamqc,
file location is configured in the parameter file):
- input_groups.txt


################################################################################
### Required external programs
### Must be executable and located within the $PATH environment.
################################################################################

Required external programs (must be executable and located within the $PATH environment):
- python3 version 3.7 or higher (executes the shell script fastq2vcf_pipeline_commands)
- R (required for the qualimap multi-bamqc)
- fastqc
- repair.sh
- fastp
- bwa-mem2 index
- bwa-mem2 mem
- samtools view
- samtools sort
- samtools index
- qualimap bamqc
- qualimap multi-bamqc


################################################################################
### How to execute the scripts
################################################################################

sbatch_fastq2vcf_pipeline_commands.sh is the master shell script which encapsulates
the specifics of the sbatch execution. This script loads the required modules
within the cheops1 environment (module load fastqc, module load bwamem2,
module load samtools) and invokes the other python shell script
fastq2vcf_pipeline_commands which does the actual work. Both scripts need
to have execution privileges, and both scripts need to be located within the
$PATH environment variable.

The master shell script is invoked as follows (for sbatch execution):
sbatch sbatch_fastq2vcf_pipeline_commands.sh [<parameter script path>]

The Python script is invoked as follows (NO sbatch execution):
fastq2vcf_pipeline_commands.sh [<parameter script path>]
