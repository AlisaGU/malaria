#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/xiao_plot/drugResistance/K13"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/xiao_plot/drugResistance/K13"
cds_each_site_codon_index="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/cds_each_site_codon_index.txt"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:

cd $working_dir
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${cds_each_site_codon_index}.sort | $bedtools intersect -a K13.bed -b stdin -wb >a
$code_dir/revert.R $working_dir/a $working_dir/b
awk '{print $4"\t"$5"\t"$6"\t"int((NR-1)/3)"\t"$7"\t"$8"\t"$9"\t"$10}' b >k13.codon.genomebase.bed
rm -rf a b

## resistance mutation
awk 'NR>3{print "PF3D7_1343700""\t"$4-1"\t"$4"\t"$0}' k13.codon.genomebase.bed | $bedtools intersect -a stdin -b mutationsRelatedResistanceBasedOnProteinIndex.bed -wa -wb | awk '{print $4"\t"$5"\t"$6"\t"$17"_"$15">"$16"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' >drugResistance.bed
