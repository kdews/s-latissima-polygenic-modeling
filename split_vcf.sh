#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=100g
#SBATCH --cpus-per-task=12
#SBATCH --time=2-0
#SBATCH -J split_vcf
#SBATCH -o %x.log

# Source conda configuration
conda_sh=~/bin/anaconda3/etc/profile.d/conda.sh
source $conda_sh

# Define input and output files
vcf=$1
sample_list=$2
[[ -f $vcf ]] && echo "Input VCF: $vcf" || { echo "No VCF provided"; exit 1; }
[[ -f $sample_list ]] && echo "Sample list: $sample_list" || \
{ echo "No list of samples provided"; exit 1; }
vcf_base=$(echo $(basename -- $vcf) | sed 's/\.vcf.*//g')
vcf_ext=$(echo $(basename -- $vcf) | sed 's/.*vcf/vcf/g')
[[ $vcf_ext == "vcf.gz" ]] && vcf_opt="--gzvcf" || vcf_opt="--vcf"
indiv=$(cat $sample_list | sed -n ${SLURM_ARRAY_TASK_ID}p)
[[ -f $remove_indvs ]] || { echo "Need list of individuals to remove"; exit 1; }
vcf_indiv=${indiv}.${vcf_base}.recode.vcf.gz

# Remove "check" individuals using VCFtools
conda activate vcftools
[[ -f $vcf_indv_filt ]] && \
echo "Using individual-filtered VCF $vcf_indv_filt" || \
{ echo "Filtering individual(s) from ${vcf}. 
Find output in $vcf_indv_filt"; \
vcftools --non-ref-ac-any 1 --indv $indiv \
$vcf_opt $vcf --recode --recode-INFO-all \
--stdout | gzip -c > $vcf_indiv; }

# hi kelly this is a test
