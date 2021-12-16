#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=100g
#SBATCH --cpus-per-task=12
#SBATCH --time=2-0
#SBATCH -J pca_s_latissima_popgen 
#SBATCH -o %x.log

# Source conda configuration
conda_sh=~/bin/anaconda3/etc/profile.d/conda.sh
source $conda_sh

# Load R module and dependencies
module load gcc/11.2.0 openblas/0.3.18 r/4.1.2
echo "Running $(Rscript --version)"

# Define input and output files
vcf=$1
[[ -f $vcf ]] && echo "Input VCF: $vcf" || { echo "No VCF provided"; exit 1; }
remove_indvs=$2
[[ -f $remove_indvs ]] || { echo "Need list of individuals to remove"; exit 1; }
vcf_base=$(echo $(basename -- $vcf) | sed 's/\.vcf.*//g')
vcf_ext=$(echo $(basename -- $vcf) | sed 's/.*vcf/vcf/g')
[[ $vcf_ext == "vcf.gz" ]] && vcf_opt="--gzvcf" || vcf_opt="--vcf"
vcf_indv_filt=${vcf_base}.indv_filt.recode.vcf.gz

# Remove "check" individuals using VCFtools
conda activate vcftools
[[ -f $vcf_indv_filt ]] && \
echo "Using individual-filtered VCF $vcf_indv_filt" || \
{ echo "Filtering individual(s) from ${vcf}. 
Find output in $vcf_indv_filt"; \
vcftools --remove $remove_indvs $vcf_opt $vcf --recode --recode-INFO-all \
--stdout | gzip -c > $vcf_indv_filt; }

# Run R to plot PCA from filtered VCF
Rscript --vanilla \
~/scripts/s-latissima-polygenic-modeling/s_latissima_popgen_pca.R \
$vcf_indv_filt

