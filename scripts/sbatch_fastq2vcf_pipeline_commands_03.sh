#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/%u/test-%j.out
#SBATCH --error=/home/group.kurse/%u/test-%j.err

# Set path for input and output directories, and reference fasta
input_dir="$HOME"/practical_3/filtered_bam
output_dir="$HOME"/practical_3/variant_calling
reference_fasta="$HOME"/practical_2/reference_fasta/Arabidopsis_thaliana.TAIR10.dna.chromosome.4_98K.fasta

module load samtools

# shellcheck disable=SC2011
samples=$(ls "$HOME"/practical_2/mapping_bwa/*Chr4_sorted.bam | xargs -n 1 basename | sed 's/\.Chr4_sorted.bam//')
for sample in $samples
do
    # Quality filtering and removing PCR duplicates with Samtools 
    # Keep only properly aligned paired-end reads (-f 3)
    # Mapping quality (-q 30)
    # Exclude secondary alignments and reads failing quality checks (-F 264)
    samtools view -b -o "$input_dir"/"$sample".Chr4_BamQualFilt-f3Q30F264.bam -f 3 -q 30 -F 264 "$HOME"/practical_2/mapping_bwa/"$sample".Chr4_sorted.bam

    # Remove duplicates (samtools rmdup)
    samtools rmdup -s "$input_dir"/"$sample".Chr4_BamQualFilt-f3Q30F264.bam "$input_dir"/"$sample".Chr4_BamQualFilt-f3Q30F264_nodup.bam
      
    # Index the filtered_nodup.bam file
    samtools index "$input_dir"/"$sample".Chr4_BamQualFilt-f3Q30F264_nodup.bam
done

module load bcftools

# List all filtered deduplicated bam files in bam_files.txt
# List all filtered deduplicated bam files in bam_files.txt
bam_files=$(find "$input_dir"/ -type f -name "*_nodup.bam")
echo "$bam_files" > "$input_dir"/bam_files.txt

# Call variants using bcftools mpileup (15 min) and call (3-5 minutes)commands
# Activate Base Alignment Quality computation (-E)
# Minimum base quality (-q 30)
bcftools mpileup -E -q 30 --threads 8 -o "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup.bcf -f "$reference_fasta" -b "$input_dir"/bam_files.txt

# Minimum calling threshold for variant alleles (-p 0.01) Variants with an allele frequency of at least 1% will be called
bcftools call -c -p 0.01 -O z --threads 8 -o "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw.vcf.gz "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup.bcf

# Basic filtering by removing all monomorphic variant sites
bcftools view -i "AC>0" "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw.vcf.gz -o "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic.vcf.gz

# Check names of samples included in the vcf file
bcftools query -l "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic.vcf.gz

# Simplify and shorten the names of samples
zcat "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic.vcf.gz | \
awk -F '\t' 'BEGIN{OFS="\t"} {if ($1 ~ /^#CHROM/) {for (i=10; i<=NF; i++) {sub(".*/", "", $i); sub("\\.Chr4_BamQualFilt-f3Q30F264_nodup\\.bam", "", $i)}} print }' | \
gzip -c > "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam.vcf.gz

# Re-check names of samples included in the vcf file
bcftools query -l "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam.vcf.gz

module unload bcftools
module unload gnu
module load vcftools

# Additional filters
vcftools --gzvcf "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam.vcf.gz --minDP 10 --minGQ 20 --minQ 30 --max-missing 0.80 --remove-indels --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz

# Count number of variants in raw vcf
zgrep -v "^#" "$output_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz | wc -l
