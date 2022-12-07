#!/bin/bash
#SBATCH -n 4
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
extract_vcf(){
    for chrom in `echo "01 02 03 04 05 06 07 08 09 10 11 12 13 14"|tr " " "\n"`
    do
    $code_dir/s4_extract_vcf.sh $chrom
    done
}

combine_vcf(){
    cd $pf6_variant_dir
    $bcftools concat -f $allvariant_list -O z -o $allvariant_file --threads 4
    $tabix $allvariant_file
    $bcftools concat -f $bisnp_list -O z -o $bisnp_file --threads 4
    $tabix $bisnp_file
    cd -
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd"
pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
allvariant_list=$pf6_variant_dir/"allvariant.list"
allvariant_file=$pf6_variant_dir/"WSEA_peak500_surr.vcf.gz"
bisnp_list=$pf6_variant_dir/"bisnp.list"
bisnp_file=$pf6_variant_dir/"WSEA_peak500_surr_bisnp.vcf.gz"

bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
tabix="/picb/evolgen/users/gushanshan/GenomeAnnotation/htslib/htslib-1.10.2/htslib_install/bin/tabix"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# extract_vcf
# combine_vcf
