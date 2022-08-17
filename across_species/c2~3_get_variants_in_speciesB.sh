#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_variants_in_speciesB() {
    chrom=$1
    maf_prefix=$2
    cd $prank_speciesB_variants_dir
    # $prank -d=$global_working_dir/${chrom}.${maf_prefix}.noFilterBlock.fasta -showanc -showevents -once -keep -o=${chrom}.anc -t=$global_working_dir/tree.nwk
    # let start_pos=$(grep -n "branch P_reichenowi" ${chrom}.anc.events | awk -F":" '{print $1}')+1
    # let end_pos=$(grep -n "branch P_bill" ${chrom}.anc.events | awk -F":" '{print $1}')-2
    # awk -v start_pos=$start_pos -v end_pos=$end_pos 'NR>=start_pos && NR<=end_pos' ${chrom}.anc.events | grep -E "insertion|deletion" >${chrom}.P_reichenowi.events
    # rm -rf ${chrom}.nonscaled_genome_pos.indel.bed
    # grep "deletion" ${chrom}.P_reichenowi.events | awk -F"." -v chrom=${chrom} '{print chrom"\t"$1-1"\t"$1"\t"0"\t"0"\t"0"\t"0"\t"1"\t"1}' >${chrom}.nonscaled_genome_pos.indel.bed
    # grep "insertion" ${chrom}.P_reichenowi.events | awk -F"." -v chrom=${chrom} '{print chrom"\t"$1-2"\t"$1-1"\t"0"\t"0"\t"0"\t"1"\t"0"\t"1}' >>${chrom}.nonscaled_genome_pos.indel.bed

    $mafFilter $global_working_dir/$chrom.${maf_prefix}.maf -minRow=2 -maxRow=2 -needComp="P_reichenowi" -reject=/dev/null | grep "P_falciparum" | awk -v chrom=$chrom '{print chrom"\t"$3"\t"$3+$4}' >${chrom}.noBill.3d7.region.bed
    # $code_dir/s2_scale_bed_and_merge_variants.r $prank_speciesB_variants_dir/${chrom}.nonscaled_genome_pos.indel.bed $global_working_dir/${chrom}.${maf_prefix}.noFilterBlock.fasta
    $bedtools subtract -a ${chrom}.scaled_genome_pos.bed -b ${chrom}.noBill.3d7.region.bed >${chrom}.scaled_genome_pos.except_no_bill_block.bed

}

export -f get_variants_in_speciesB
# VARIABLE NAMING TODO:
export code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/across_species"
export global_working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf"
export prank_speciesB_variants_dir=$global_working_dir/prank_speciesB_variants
maf_prefix="P_falciparum.P_billcollinsi.P_reichenowi"

parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export prank="/picb/evolgen/users/gushanshan/software/PRANK/bin/prank"
export mafFilter="/picb/evolgen/users/gushanshan/software/ucsc_pairwiseGenomeAlignmrent/data/bin/mafFilter"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
mkdir -p $prank_speciesB_variants_dir
$parallel -j 14 get_variants_in_speciesB ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix
