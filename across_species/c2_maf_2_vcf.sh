#!/bin/bash
#SBATCH -n 60
#SBATCH -p hpc
#SBATCH -o %x-%j.out
#SBATCH --mem=300G
# FUNCTIONS TODO:
maf2fasta_species3() {
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
    # maf_dir=$working_dir/P_falciparum.P_reichenowi
    # mkdir -p $maf_dir
    # cd $maf_dir
    # # $mafFilter ../$chrom.${maf_prefix}.maf -minRow=2 -needComp="P_reichenowi" -reject=${chrom}.reject2.txt | sed '/P_praefalciparum/d' >$chrom.P_falciparum.P_reichenowi.maf
    # $mafFilter ../$chrom.${maf_prefix}.maf -minRow=2 -needComp="P_reichenowi" -reject=${chrom}.reject2.txt | sed '/P_billcollinsi/d' >$chrom.P_falciparum.P_reichenowi.maf

    # $code_dir/s2_maf2fasta_species2.R "P_falciparum.P_reichenowi" $chrom.P_falciparum.P_reichenowi.maf

    # cat $chrom.P_falciparum.fasta $chrom.P_reichenowi.fasta >$chrom.P_falciparum.P_reichenowi.fasta
    # rm -rf $chrom.P_falciparum.fasta $chrom.P_reichenowi.fasta

    # 寻找包含两个物种（P_falciparum、P_praefalciparum）的block
    # maf_dir=$working_dir/P_falciparum.P_praefalciparum
    # mkdir -p $maf_dir
    # cd $maf_dir
    # $mafFilter ../$chrom.${maf_prefix}.maf -minRow=2 -needComp="P_praefalciparum" -reject=${chrom}.reject2.txt | sed '/P_reichenowi/d' >$chrom.P_falciparum.P_praefalciparum.maf
    # $code_dir/s2_maf2fasta_species2.R "P_falciparum.P_praefalciparum" $chrom.P_falciparum.P_praefalciparum.maf
    # cat $chrom.P_falciparum.fasta $chrom.P_praefalciparum.fasta >$chrom.P_falciparum.P_praefalciparum.fasta
    # rm -rf $chrom.P_falciparum.fasta $chrom.P_praefalciparum.fasta
}

fasta2vcf_species3() {
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
        # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==$11 && $10=="0" &&$10!=$12 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species3.raw.vcf
        awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=3{print};NR==4{$11=null;print};NR>4 && $10==$11 && $10!=$12 {$11=null;$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species3.raw.include_indel.vcf

    fi

    # $code_dir/s2_scale_vcf_ordinate.r ${chrom}.species3.raw.vcf ${chrom}.${maf_prefix}.species3.fasta
    # rm -rf ${chrom}.species3.raw.vcf ${chrom}.vcf.intermediate

    $code_dir/s2_scale_vcf_ordinate_include_insertion.r ${chrom}.species3.raw.include_indel.vcf $chrom.*.species3.fasta

    rm -rf ${chrom}.species3.raw.include_indel.vcf ${chrom}.vcf.intermediate

    # 两个物种（P_falciparum、P_reichenowi）的序列
    # maf_dir=$working_dir/P_falciparum.P_reichenowi
    # cd $maf_dir
    # $snp_sites ${chrom}.P_falciparum.P_reichenowi.fasta -b -v >${chrom}.vcf.intermediate
    # last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    # non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    # if [ $non_comment_line_count -eq $last_pos ]; then
    #     # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==0 && $11==1 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species2.raw.vcf
    #     awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10!="." && $10!=$11 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species2.raw.include_indel.vcf

    # fi

    # # $code_dir/s2_scale_vcf_ordinate.r P_falciparum.P_reichenowi/${chrom}.species2.raw.vcf P_falciparum.P_reichenowi/${chrom}.P_falciparum.P_reichenowi.fasta
    # $code_dir/s2_scale_vcf_ordinate_include_insertion.r P_falciparum.P_reichenowi/${chrom}.species2.raw.include_indel.vcf P_falciparum.P_reichenowi/${chrom}.P_falciparum.P_reichenowi.fasta

    # rm -rf ${chrom}.species2.raw.include_indel.vcf ${chrom}.vcf.intermediate

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

maf2fasta_species3_noFilterBlock() {
    chrom=$1
    maf_prefix=$2
    $mafFilter ${chrom}.${maf_prefix}.maf -needComp="P_falciparum" -minRow=1 -reject=${chrom}.reject.noP_falciparum.txt >${chrom}.${maf_prefix}.includingP_falciparum.maf

    value="1"
    while read line; do
        line_number=$(grep -n "$line" ${chrom}.${maf_prefix}.maf | awk -F":" '{print $1}')
        let single_block_start=line_number-1
        let former_last_block_end=line_number-3
        let latter_first_block_start=line_number+2
        value+=" $former_last_block_end $single_block_start $line_number $latter_first_block_start"
    done < <(grep "^s " ${chrom}.reject.noP_falciparum.txt)
    value+=" $(wc -l ${chrom}.${maf_prefix}.maf | awk '{print $1}')"

    let file_count=$(echo $value | tr " " "\n" | wc -l)/2
    for file_order in $(seq 1 $file_count); do
        file_start=$(echo $value | tr " " "\n" | awk -v file_order=$file_order 'NR==2*file_order-1')
        file_end=$(echo $value | tr " " "\n" | awk -v file_order=$file_order 'NR==2*file_order')
        awk -v file_start=$file_start -v file_end=$file_end 'NR>=file_start && NR<=file_end' ${chrom}.${maf_prefix}.maf >${chrom}.file${file_order}

        let oushu=$file_order%2
        if [ $oushu -eq 0 ]; then
            awk '$1~/^s/ {$5="+";print}' ${chrom}.file${file_order} >${chrom}.file${file_order}.1
            mv ${chrom}.file${file_order}.1 ${chrom}.file${file_order}
        fi
    done

    $code_dir/s2_maf2fasta_noFilterBlock.R $chrom
    cat ${chrom}.P_falciparum.noFilterBlock.fasta ${chrom}.P_billcollinsi.noFilterBlock.fasta ${chrom}.P_reichenowi.noFilterBlock.fasta >${chrom}.P_falciparum.P_billcollinsi.P_reichenowi.noFilterBlock.fasta
    rm -rf ${chrom}.P_falciparum.noFilterBlock.fasta ${chrom}.P_billcollinsi.noFilterBlock.fasta ${chrom}.P_reichenowi.noFilterBlock.fasta ${chrom}.file* ${chrom}.${maf_prefix}.includingP_falciparum.maf
}

fasta2vcf_species3_noFilterBlock_bedLatter() {
    chrom=$1
    maf_prefix=$2

    $snp_sites ${chrom}.${maf_prefix}.noFilterBlock.fasta -b -v >${chrom}.vcf.intermediate
    last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    if [ $non_comment_line_count -eq $last_pos ]; then
        # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==$12 && $10=="0" &&$10!=$11 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species3.raw.vcf
        # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10==$11 && $10=="0" &&$10!=$12 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species3.raw.vcf
        awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=3{print};NR==4{$11=null;print};NR>4 && $10==$11 && $10!=$12 {$11=null;$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.noFilterBlock.raw.include_indel.vcf

    fi

    $code_dir/s2_scale_vcf_ordinate_include_insertion.r ${chrom}.noFilterBlock.raw.include_indel.vcf $chrom.*.noFilterBlock.fasta

    rm -rf ${chrom}.noFilterBlock.raw.include_indel.vcf ${chrom}.vcf.intermediate
}

fasta2vcf_species3_noFilterBlock_bedFirst() {
    chrom=$1
    maf_prefix=$2
    cd $working_dir
    # $snp_sites ${chrom}.${maf_prefix}.noFilterBlock.fasta -b -v >${chrom}.vcf.intermediate

    # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10!="." && $11!="." && $12!="." {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.vcf.intermediate1

    # $code_dir/s2_get_mutation_load_for_each_chrom_each_site.R $working_dir/${chrom}.vcf.intermediate1
    # $code_dir/s2_get_mutation_load_for_each_chrom_each_site_Pb_as_ref.R $working_dir/${chrom}.vcf.intermediate1
    # $code_dir/s2_get_sutiable_site_for_indel.R $working_dir/${chrom}.vcf.intermediate1
    # $code_dir/s2_scale_bed.r $chrom.vcf.indel_possible.bed $chrom.$maf_prefix.noFilterBlock.fasta
    $code_dir/s2_scale_bed.r $chrom.intermediate1_only_indel_polysite.bed $chrom.$maf_prefix.noFilterBlock.fasta
    mv $chrom.intermediate1_only_indel_polysite.scaled_genome_pos.bed $chrom.twoOutgroup_nogap_only_indel.bed
    # awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR>4 && $10==$11 && $10!=$12 {$11=null;$1=chrom;print $1"\t"$2-1"\t"$2}' ${chrom}.vcf.intermediate | $bedtools intersect -a ${chrom}.intermediate1_polysite.bed -b stdin >${chrom}.noFilterBlock.raw.include_indel.bed

    # mv ${chrom}.intermediate1_polysite.bed ${chrom}.allvariant.noFilterBlock.raw.include_indel.bed

    # $code_dir/s2_scale_bed_ordinate_include_insertion.r ${chrom}.allvariant.noFilterBlock.raw.include_indel.bed $chrom.$maf_prefix.noFilterBlock.fasta
    # $code_dir/s2_scale_bed_ordinate_include_insertion.r ${chrom}.noFilterBlock.raw.include_indel.bed $chrom.$maf_prefix.noFilterBlock.fasta
    # $bedtools intersect -a ${chrom}.allvariant.noFilterBlock.include_indel.bedFirst.bed -b ${chrom}.vcf.indel_possible.scaled_genome_pos.bed >${chrom}.onlyIndel.noFilterBlock.include_indel.bedFirst.bed

    # rm -rf ${chrom}.vcf.intermediate ${chrom}.vcf.intermediate1 $chrom.vcf.indel_possible.bed ${chrom}.vcf.indel_possible.scaled_genome_pos.bed ${chrom}.noFilterBlock.raw.include_indel.bed ${chrom}.allvariant.noFilterBlock.raw.include_indel.bed
}

maf2vcf_species3_noFilterBlock_bedFirst_Li() {
    chrom=$1
    maf_prefix=$2
    cd $working_dir

    # $mafFilter ${chrom}.${maf_prefix}.maf -needComp="P_billcollinsi" -minRow=2 -maxRow=3 -reject=${chrom}.no.P_billcollinsi.maf >${chrom}.${maf_prefix}.include.P_billcollinsi.maf
    # $code_dir/s2_maf2fasta.R ${chrom}.${maf_prefix}.include.P_billcollinsi.maf

    # cmd="cat"
    # for i in $(seq 1 3); do
    #     cmd+=" $chrom.$(echo $maf_prefix | awk -F "\." -v i=$i '{print $i}' 2>/dev/null).fasta"
    # done
    # cmd+=" >$chrom.include_P_billcollinsi.fasta"
    # eval $cmd

    # cmd="rm -rf "
    # for i in $(seq 1 3); do
    #     cmd+=" $chrom.$(echo $maf_prefix | awk -F "\." -v i=$i '{print $i}' 2>/dev/null).fasta"
    # done
    # eval $cmd

    # $snp_sites ${chrom}.include_P_billcollinsi.fasta -b -v | awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4{$1=chrom;print $0}' >${chrom}.vcf.include_P_billcollinsi.intermediate

    ##从${chrom}.vcf.include_P_billcollinsi.intermediate中挑选变异，有两种方式，均由s2_get_mutation_load_for_each_chrom_each_site_Pb_as_ref_Li.R文件产生：
    ##（1）考虑标准株是gap的情况，这时，indel位置不准确，结果文件是${chrom}.include_P_billcollinsi.intermediate_only_indel_polysite.bed
    ##（2）不考虑标准株是gap，定位准确，               结果文件是${chrom}.include_P_billcollinsi.intermediate_only_indel_polysite_standardStrainBaseExist.bed
    $code_dir/s2_get_mutation_load_for_each_chrom_each_site_Pb_as_ref_Li.R $working_dir/${chrom}.vcf.include_P_billcollinsi.intermediate
    # $code_dir/s2_scale_bed.r $chrom.include_P_billcollinsi.intermediate_only_indel_polysite.bed $chrom.include_P_billcollinsi.fasta
    $code_dir/s2_scale_bed.r ${chrom}.include_P_billcollinsi.intermediate_only_indel_polysite_standardStrainBaseExist.bed $chrom.include_P_billcollinsi.fasta

    ##从block角度出发，选择有外群的block
    # awk -v chrom=$chrom '$2~/P_falciparum/{print chrom"\t"$3"\t"$3+$4}' $chrom.P_falciparum.P_billcollinsi.P_reichenowi.include.P_billcollinsi.maf >${chrom}.include_P_billcollinsi.bed
    ##相比于block角度，进一步深化，挑选外群出现的位置。舍弃block中有一部分没有外群的地方
    # $code_dir/s2_get_outgroupA_region.r $working_dir/${chrom}.include_P_billcollinsi.fasta $chrom # 搜索P_billcollinsi出现的区域，不考虑其他物种
}

split_by_chromosome() {
    ChrName=$1
    filename=$2
    chrAlignStart=$(grep -n "s P_falciparum.${ChrName} " $filename | head -1 | cut -d ":" -f 1)
    chrAlignStart=$(expr $chrAlignStart - 1)
    chrAlignEnd=$(grep -n "s P_falciparum.${ChrName} " $filename | tail -1 | cut -d ":" -f 1)
    chrAlignEnd=$(expr $chrAlignEnd + 2)
    sed -n "${chrAlignStart},${chrAlignEnd}p" $filename >./${ChrName}.$(basename $filename)
}

maf2fasta_species2() {
    chrom=$1
    cd $maf_dir
    $code_dir/s2_maf2fasta_species2.R "P_falciparum.P_reichenowi" $chrom.P_falciparum.P_reichenowi.maf

    cat $chrom.P_falciparum.fasta $chrom.P_reichenowi.fasta >$chrom.P_falciparum.P_reichenowi.fasta
    rm -rf $chrom.P_falciparum.fasta $chrom.P_reichenowi.fasta
}

fasta2vcf_species2() {
    chrom=$1

    # 两个物种（P_falciparum、P_reichenowi）的序列
    cd $maf_dir
    $snp_sites ${chrom}.P_falciparum.P_reichenowi.fasta -b -v >${chrom}.vcf.intermediate
    last_pos=$(tail -n 1 ${chrom}.vcf.intermediate | awk '{print $2}')
    non_comment_line_count=$(grep -v "^#" ${chrom}.vcf.intermediate | wc -l)

    if [ $non_comment_line_count -eq $last_pos ]; then
        awk -v chrom=$chrom 'BEGIN{OFS="\t"};NR<=4{print};NR>4 && $10!="." && $10!=$11 {$1=chrom;print}' ${chrom}.vcf.intermediate >${chrom}.species2.raw.include_indel.vcf
    fi

    $code_dir/s2_scale_vcf_ordinate_include_insertion.r P_falciparum.P_reichenowi/${chrom}.species2.raw.include_indel.vcf P_falciparum.P_reichenowi/${chrom}.P_falciparum.P_reichenowi.fasta

    rm -rf ${chrom}.species2.raw.include_indel.vcf ${chrom}.vcf.intermediate
}
export -f maf2fasta_species2 split_by_chromosome maf2fasta_species3_noFilterBlock fasta2vcf_species3_noFilterBlock_bedFirst maf2vcf_species3_noFilterBlock_bedFirst_Li

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
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
cd $working_dir

# $maffilter param=option

# from 3-ways alignment
# $parallel maf2fasta_species3 ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix

# $parallel fasta2vcf_species3 ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix

# $parallel maf2fasta_species3_noFilterBlock ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix
# $parallel fasta2vcf_species3_noFilterBlock_bedLatter ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix
# $parallel -j 14 fasta2vcf_species3_noFilterBlock_bedFirst ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix
$parallel -j 14 maf2vcf_species3_noFilterBlock_bedFirst_Li ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: $maf_prefix

# from 2-ways alignment
# export maf_dir=$working_dir/P_falciparum.P_reichenowi
# mkdir -p $maf_dir
# cd $maf_dir
# $parallel -j 4 split_by_chromosome ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3 ::: P_falciparum.P_reichenowi.maf
# $parallel -j 4 maf2fasta_species2 ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3
# $parallel -j 4 fasta2vcf_species2 ::: Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3
