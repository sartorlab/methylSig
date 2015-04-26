# All paths are relative to working environment on our lab's compute cluster
cd ~/latte/vcf/

# Shell script to pull out C > T SNPs from ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz

# Filter the sites.vcf so that AF[0]>0.05
# NOTE: If AF > 0.5, then the 'alternative' allele is actually the main one
~/latte/apps/bcftools-1.2/bin/bcftools filter -i 'AF[0]>0.05' --output-type=v --output=filtered_AF_0.05.vcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf

# Grab only C > T (forward strand) and G > A (reverse strand)
grep -P '(C\tT)' filtered_AF_0.05.vcf > filtered_AF_0.05_CT_forward.vcf
grep -P '(G\tA)' filtered_AF_0.05.vcf > filtered_AF_0.05_CT_reverse.vcf

# Concatenate the two files
cat filtered_AF_0.05_CT_forward.vcf filtered_AF_0.05_CT_reverse.vcf > filtered_AF_0.05_CT.vcf

# Grab only the SNP variants
grep -P '(VT=SNP)' filtered_AF_0.05_CT.vcf > filtered_AF_0.05_CT_SNPs.vcf

# NOTE VCFs are 1-based
