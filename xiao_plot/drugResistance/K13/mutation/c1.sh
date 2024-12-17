#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
vcf_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/pf6"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome
faToTwoBit PlasmoDB-36_Pfalciparum3D7_Genome.fasta PlasmoDB-36_Pfalciparum3D7_Genome.2bit
less -S PlasmoDB-36_Pfalciparum3D7_Genome.fasta | grep ">" | awk '{print $1}' | sed 's/^>//g' >pf36.chromName.txt

cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13/variant_annotation
zcat $vcf_dir/Pf_60_public_Pf3D7_13_v3.final.vcf.gz | head -n 500 | grep "#CHROM" >header
zcat $vcf_dir/Pf_60_public_Pf3D7_13_v3.final.vcf.gz | $bedtools intersect -a stdin -b K13.bed >K13.vcf
cat header K13.vcf >K13.complete.vcf
