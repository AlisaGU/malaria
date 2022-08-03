#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
maf2fasta() {
    chrom=$1
    maf_prefix=$2
    # 寻找包含P_falciparum的block
    # sed -i '1i##maf version=1' $chrom.${maf_prefix}.maf
    # $mafFilter $chrom.${maf_prefix}.maf -speciesFilter=$working_dir/species.lst -minRow=1 -reject=${chrom}.reject.txt >$chrom.${maf_prefix}.filter.maf
    # $code_dir/s2_maf2fasta.R $chrom.${maf_prefix}.filter.maf
    # cat $chrom.P_falciparum.fasta $chrom.P_praefalciparum.fasta $chrom.P_reichenowi.fasta >$chrom.${maf_prefix}.fasta
    # rm -rf $chrom.P_falciparum.fasta $chrom.P_praefalciparum.fasta $chrom.P_reichenowi.fasta

    # 寻找三个物种都有的block
    sed -i '1i##maf version=1' $chrom.${maf_prefix}.maf
    $mafFilter $chrom.${maf_prefix}.maf -minRow=3 -maxRow=3 -reject=${chrom}.reject3.txt >$chrom.${maf_prefix}.filter3.maf
    $code_dir/s2_maf2fasta.R $chrom.${maf_prefix}.filter3.maf

    cmd="cat"
    for i in $(seq 1 3); do
        cmd+=" $chrom.$(echo $maf_prefix | awk -F "\." -v i=$i '{print $i}' 2>/dev/null).fasta"
    done
    cmd+=" >$chrom.${maf_prefix}.species3.fasta"
    eval $cmd

    cmd="rm -rf "
    for i in $(seq 1 3); do
        cmd+=" $chrom.$(echo $maf_prefix | awk -F "\." -v i=$i '{print $i}' 2>/dev/null).fasta"
    done
    eval $cmd

    # 寻找包含两个物种（P_falciparum、P_reichenowi）的block
    maf_dir=$working_dir/P_falciparum.P_reichenowi
    mkdir -p $maf_dir
    cd $maf_dir
    # $mafFilter ../$chrom.${maf_prefix}.maf -minRow=2 -needComp="P_reichenowi" -reject=${chrom}.reject2.txt | sed '/P_praefalciparum/d' >$chrom.P_falciparum.P_reichenowi.maf
    $mafFilter ../$chrom.${maf_prefix}.maf -minRow=2 -needComp="P_reichenowi" -reject=${chrom}.reject2.txt | sed '/P_billcollinsi/d' >$chrom.P_falciparum.P_reichenowi.maf

    $code_dir/s2_maf2fasta_species2.R "P_falciparum.P_reichenowi" $chrom.P_falciparum.P_reichenowi.maf

    cat $chrom.P_falciparum.fasta $chrom.P_reichenowi.fasta >$chrom.P_falciparum.P_reichenowi.fasta
    rm -rf $chrom.P_falciparum.fasta $chrom.P_reichenowi.fasta

    # 寻找包含两个物种（P_falciparum、P_praefalciparum）的block
    # maf_dir=$working_dir/P_falciparum.P_praefalciparum
    # mkdir -p $maf_dir
    # cd $maf_dir
    # $mafFilter ../$chrom.${maf_prefix}.maf -minRow=2 -needComp="P_praefalciparum" -reject=${chrom}.reject2.txt | sed '/P_reichenowi/d' >$chrom.P_falciparum.P_praefalciparum.maf
    # $code_dir/s2_maf2fasta_species2.R "P_falciparum.P_praefalciparum" $chrom.P_falciparum.P_praefalciparum.maf
    # cat $chrom.P_falciparum.fasta $chrom.P_praefalciparum.fasta >$chrom.P_falciparum.P_praefalciparum.fasta
    # rm -rf $chrom.P_falciparum.fasta $chrom.P_praefalciparum.fasta
}

fasta2vcf() {
    chrom=$1
    maf_prefix=$2

    #只要P_falciparum有序列就保留
    # $snp_sites ${chrom}.${maf_prefix}.fasta -b -v >${chrom}.vcf.intermediate
    # last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    # non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    # if [ $non_comment_line_count -eq $last_pos ]; then
    #     awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10!="." && $11!="."&& $12!="."{$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.raw.vcf
    # fi

    # # java -jar /picb/evolgen/users/gushanshan/software/jvarkit/jvarkit/dist/msa2vcf.jar ${chrom}.${maf_prefix}.fasta -R $chrom >${chrom}.raw.vcf
    # $code_dir/s2_scale_vcf_ordinate.r ${chrom}.raw.vcf ${chrom}.${maf_prefix}.fasta
    # rm -rf ${chrom}.raw.vcf ${chrom}.vcf.intermediate

    # 3个物种都得有
    $snp_sites ${chrom}.${maf_prefix}.species3.fasta -b -v >${chrom}.vcf.intermediate
    last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    if [ $non_comment_line_count -eq $last_pos ]; then
        # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==$12 && $10=="0" &&$10!=$11 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species3.raw.vcf
        awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==$11 && $10=="0" &&$10!=$12 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species3.raw.vcf
    fi

    $code_dir/s2_scale_vcf_ordinate.r ${chrom}.species3.raw.vcf ${chrom}.${maf_prefix}.species3.fasta
    rm -rf ${chrom}.species3.raw.vcf ${chrom}.vcf.intermediate

    # 两个物种（P_falciparum、P_reichenowi）的序列
    maf_dir=$working_dir/P_falciparum.P_reichenowi
    cd $maf_dir
    $snp_sites ${chrom}.P_falciparum.P_reichenowi.fasta -b -v >${chrom}.vcf.intermediate
    last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    if [ $non_comment_line_count -eq $last_pos ]; then
        awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==0 && $11==1 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species2.raw.vcf
    fi

    $code_dir/s2_scale_vcf_ordinate.r P_falciparum.P_reichenowi/${chrom}.species2.raw.vcf P_falciparum.P_reichenowi/${chrom}.P_falciparum.P_reichenowi.fasta
    rm -rf ${chrom}.species2.raw.vcf ${chrom}.vcf.intermediate

    # 两个物种（P_falciparum、P_praefalciparum）的序列
    # maf_dir=$working_dir/P_falciparum.P_praefalciparum
    # cd $maf_dir
    # $snp_sites ${chrom}.P_falciparum.P_praefalciparum.fasta -b -v >${chrom}.vcf.intermediate
    # last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    # non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    # if [ $non_comment_line_count -eq $last_pos ]; then
    #     awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==0 && $11==1 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species2.raw.vcf
    # fi

    # $code_dir/s2_scale_vcf_ordinate.r P_falciparum.P_praefalciparum/${chrom}.species2.raw.vcf P_falciparum.P_praefalciparum/${chrom}.P_falciparum.P_praefalciparum.fasta
    # rm -rf ${chrom}.species2.raw.vcf ${chrom}.vcf.intermediate
}

export -f maf2fasta fasta2vcf

# VARIABLE NAMING TODO:
export code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/across_species"
export working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/maf2vcf"
# maf_prefix="P_falciparum.P_reichenowi.P_praefalciparum"
maf_prefix="P_falciparum.P_billcollinsi.P_reichenowi"

maffilter="/picb/evolgen/users/gushanshan/software/mafFilter/MafFilter/maffilter"
export mafFilter="/picb/evolgen/users/gushanshan/software/ucsc_pairwiseGenomeAlignmrent/data/bin/mafFilter"
msa2vcf="/picb/evolgen/users/gushanshan/software/jvarkit/jvarkit/dist/msa2vcf.jar"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export snp_sites="/picb/evolgen/users/gushanshan/software/snp-sites/snp_sites_install/bin/snp-sites"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $working_dir

# $maffilter param=option

$parallel maf2fasta ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix

$parallel fasta2vcf ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix
