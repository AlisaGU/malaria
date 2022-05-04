#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --mem=15G

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
chrom=chrom_template
export pf6_variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
export bcftools="/picb/evolgen/users/gushanshan/GenomeAnnotation/bcftools/bcftools-1.10.2/bcftools_install/bin/bcftools"
distant_african_strain="PF0870-C" # 非洲株里距离3D7最远的菌株
reference_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7"
anno_dir=${reference_dir}/anno
two_outgroup_working_dir=$anno_dir/3D7_distant_africa_strain_consistent_region_in_core
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
cd $two_outgroup_working_dir
variant_file=$pf6_variant_dir/"Pf_60_public_Pf3D7_${chrom}_v3.final.vcf.gz"
$bcftools view --include 'FILTER="PASS"' -s $distant_african_strain $variant_file | grep -v "^#" | awk 'BEGIN{OFS="\t"}{for(i=1;i<=2;i++){printf("%s%s",$i,OFS)}for(i=10;i<=10;i++){split($i,a,":");split(a[1],b,/[|/]/);printf("%s",b[2])};printf(ORS)}' | awk '$3!=0{print $1"\t"$2-1"\t"$2}' >${chrom}.inconsistent_region_in_core.bed
