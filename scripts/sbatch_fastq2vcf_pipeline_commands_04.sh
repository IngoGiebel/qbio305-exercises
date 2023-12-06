#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/%u/%j.out
#SBATCH --error=/home/group.kurse/%u/%j.err


# Load required modules -------------------------------------------------------

module load bcftools


# Environment vars for the data directories -----------------------------------

reference_fasta_file="$HOME"/practical_2/reference_fasta/Arabidopsis_thaliana.TAIR10.dna.chromosome.4_98K.fasta
filtered_bam_dir="$HOME"/practical_3/filtered_bam
variant_calling_dir="$HOME"/practical_3/variant_calling


# Output some information about the VCF file ----------------------------------

echo Sample list:
bcftools query -l "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz

echo Number of samples:
bcftools query -l "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz | wc -l

echo List of positions:
bcftools query -f "%POS\n" "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz | head -n 10

echo List of positions and alleles:
bcftools query -f "%CHROM %POS %REF %ALT\n" "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz | head -n 10

echo Per-sample tags:
bcftools query -f "%CHROM %POS[\t%GT\t%PL]\n" "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz | head -n 10


# Create VCF file with DP /AD in FORMAT colum ---------------------------------

# Call variants using bcftools mpileup (15 min) and call (3-5 minutes)commands
# Activate Base Alignment Quality computation (-E)
# Minimum base quality (-q 30)
bcftools mpileup -E -q 30 --threads 8 -o "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP.bcf --annotate FORMAT/AD,FORMAT/DP -f "$reference_fasta_file" -b "$filtered_bam_dir"/bam_files.txt

# Minimum calling threshold for variant alleles (-p 0.01) Variants with an allele frequency of at least 1% will be called
bcftools call -c -p 0.01 -O z --threads 8 -o "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_raw.vcf.gz "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP.bcf

# Basic filtering by removing all monomorphic variant sites
bcftools view -i "AC>0" "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_raw.vcf.gz -o "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_OnlyPolymorphic.vcf.gz

# Simplify and shorten the names of samples
zcat "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_OnlyPolymorphic.vcf.gz | \
awk -F '\t' 'BEGIN{OFS="\t"} {if ($1 ~ /^#CHROM/) {for (i=10; i<=NF; i++) {sub(".*/", "", $i); sub("\\.Chr4_BamQualFilt-f3Q30F264_nodup\\.bam", "", $i)}} print }' | \
gzip -c > "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_OnlyPolymorphic_shortNam.vcf.gz

module unload bcftools
module unload gnu
module load vcftools

# Additional filters
vcftools --gzvcf "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_OnlyPolymorphic_shortNam.vcf.gz --minDP 10 --minGQ 20 --minQ 30 --max-missing 0.80 --remove-indels --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > "$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz
