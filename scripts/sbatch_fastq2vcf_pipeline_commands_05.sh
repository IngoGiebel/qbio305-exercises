#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/%u/%j.out
#SBATCH --error=/home/group.kurse/%u/%j.err


# Environment vars for the data file and directories --------------------------

reference_fasta_file="$HOME"/practical_2/reference_fasta/Arabidopsis_thaliana.TAIR10.dna.chromosome.4_98K.fasta
filtered_bam_dir="$HOME"/practical_3/filtered_bam
variant_calling_dir="$HOME"/practical_3/variant_calling
vcf_file="$variant_calling_dir"/all_samples_Chr4_BamQualFilt-f3Q30F264_nodup_DP_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz
plink_dir="$HOME"/practical_5/plink
admix_dir="$HOME"/practical_5/admixture
admix_file=a.thaliana_admix
stacks_dir="$HOME"/practical_5/stacks


# Load required modules -------------------------------------------------------

module load stacks/2.65


# PCA -------------------------------------------------------------------------

# Change working directory to the plink subdirectory
cd $plink_dir

# Perform linkage pruning - i.e. identify prune sites
plink --vcf $vcf_file --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out a_thaliana

# Prune and create PCA
plink --vcf $vcf_file --double-id --allow-extra-chr --set-missing-var-ids @:# --extract a_thaliana.prune.in --make-bed --pca --out "$plink_dir"/a_thaliana


# Admixture analysis ----------------------------------------------------------

# Change working directory to the admixture subdirectory
cd $admix_dir

# Generate the input file in plink format
plink --double-id --vcf $vcf_file --make-bed --out $admix_file --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes.
# We will thus just exchange the first column by 0
awk '{$1="0";print $0}' "$admix_file".bim > "$admix_file".bim.tmp
mv "$admix_file".bim.tmp "$admix_file".bim

# Let’s now run it in a for loop with K=2 to K=10 and direct the output into log files
for i in {1..10}
do
  admixture --cv "$admix_file".bed $i > log${i}.out
done

# To identify the best value of k clusters which is the value with lowest cross-validation error, we need to collect the cv errors. 
# Below are three different ways to extract the number of K and the CV error for each corresponding K. 
# Like we said at the start of the course, there are many ways to achieve the same thing in bioinformatics!
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > "$admix_file".cv.error
cat "$admix_file".cv.error

# 1:0.22746
# 1 0.22574
# 2 0.21294
# 3 0.20820
# 4 0.20166
# 5 0.20605
# 6 0.20924
# 7 0.21264
# 8 0.21588
# 9 0.22131

# It seems cluster 4 is the best!

# Find the optimal number of clusters
# Collect the cross validation information obtained from the all the log files.
grep -h CV log*.out>cross_validation.txt


# Prepare STRUCTURE input file ------------------------------------------------

# Estimate population specific stats: Expected (He) and Observed heterozygosity (Ho), Pi, FIS
populations -V $vcf_file -M "$stacks_dir"/pop.txt -t 10 --structure -O "$stacks_dir"
