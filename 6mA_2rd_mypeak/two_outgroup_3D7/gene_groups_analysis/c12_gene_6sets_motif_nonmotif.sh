#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
locate_motif_bed_in_gene() {
    cd $working_dir
    for gene_set in $(echo $gene_sets | tr " " "\n"); do
        cd $gene_set
        ## 取出基因对应的序列
        # $seqkit subseq --bed ${gene_set}_genes.bed $genome -o ${gene_set}.fa
        ## 计算基因的长度
        # grep "^>" ${gene_set}.fa | sed 's/>//g' | awk '{split($1,a,"_|-|:");s=0;e=a[5]-a[4]+1;print $1"\t"s"\t"e}' >${gene_set}.gene.len
        ## 找到基因上可能存在的motif的位置
        cat ${gene_set}.fa | $seqkit locate -i -d -f $motif_pattern --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e}' >${gene_set}.motif.bed
        cd -
    done
}

get_mutation_count_for_motif_nonmotif_each_gene_set() {
    gene_set=$1
    var_type=$2

    mutation_density_dir=$mutation_density_global_dir/${gene_set}
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count
    cd $working_dir/$gene_set
    cat ${gene_set}_genes.bed | while read line; do
        motif_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.motif.bed | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $motif_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA" "NA" "NA")"
        elif [ $motif_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            ##motif region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+="$(echo $gene_name "$result")"
            ##nonmotif region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b ${gene_set}.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
        fi
        echo $value >>$mean_variant_count
    done
}

get_mutation_count_for_Mergedmotif_nonmotif_each_gene_set() {
    gene_set=$1
    var_type=$2

    mutation_density_dir=$mutation_density_global_dir/${gene_set}
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_Merged
    rm -rf $mean_variant_count
    cd $working_dir/$gene_set
    cat ${gene_set}_genes.bed | while read line; do
        motif_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.merged.motif.bed | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $motif_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA" "NA" "NA")"
        elif [ $motif_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            ##motif region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.merged.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+="$(echo $gene_name "$result")"
            ##nonmotif region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b ${gene_set}.merged.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
        fi
        echo $value >>$mean_variant_count
    done
}

get_mutation_count_for_peakmotif_nonpeakmotif_other_each_gene_set() {
    gene_set=$1
    var_type=$2

    mutation_density_dir=$mutation_density_global_dir/${gene_set}
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count_peakmotif_nonpeakmotif_other
    rm -rf $mean_variant_count
    cd $working_dir/$gene_set
    cat ${gene_set}_genes.bed | while read line; do
        peakmotif_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.motif.bed | $bedtools intersect -a stdin -b $peak_loc | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $peakmotif_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA")"
        elif [ $peakmotif_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            ##peak motif region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | $bedtools intersect -a stdin -b ${gene_set}.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+="$(echo $gene_name "$result")"
            ##peak nonmotif region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | $bedtools subtract -a stdin -b ${gene_set}.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+="$(echo " $result")"
            ##control motif region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_loc | $bedtools intersect -a stdin -b ${gene_set}.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
            ##control nonmotif region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_loc | $bedtools subtract -a stdin -b ${gene_set}.motif.bed | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
        fi
        echo $value >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region() {
    chrom=$1
    var_type=$2
    region_bed=$3

    chrom_number=$(echo $chrom | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(echo $"$region_bed" | awk '{print $3-$2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        # In this case, the ":" is just a "do-nothing" place holder
        :
    else
        chrom_mutation_count=0
    fi

    echo "$chrom_mutation_count $chrom_region_len"
}

summary_gene_ATCG_content() {
    gene_set=$1
    cd $working_dir/$gene_set
    $seqkit fx2tab -B "A" -B "T" -B "C" -B "G" -n ${gene_set}.fa | awk -F" " '{print $2"\t"$3"\t"$4"\t"$5"\t"$6}' >${gene_set}.base_content
}

get_merged_motif_bases_and_motif_count() {
    gene_set=$1
    cd $working_dir/$gene_set
    gene_motif_info=${gene_set}_motif_info
    rm -rf $gene_motif_info

    cat ${gene_set}_genes.bed | while read line; do
        motif_base_count=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.merged.motif.bed | awk '{print $3-$2}' | awk '{s+=$1};END{print s}')
        if [ $motif_base_count -gt 0 ] 2>/dev/null; then
            :
        else
            motif_base_count=0
        fi
        motif_count=$(echo $"$line" | $bedtools intersect -a stdin -b ${gene_set}.motif.bed | wc -l)
        gene_name=$(echo $line | awk '{print $4}')
        gene_len=$(echo $line | awk '{print $3-$2}')
        echo $gene_name" "$gene_len" "$motif_count" "$motif_base_count >>$gene_motif_info
    done
}

explore_why() {
    gene_set=$1
    cd $working_dir/$gene_set
    sort -k1,1 -k2,2n ${gene_set}.motif.bed | $bedtools merge -i stdin >${gene_set}.merged.motif.bed
}
export -f get_mutation_count_for_motif_nonmotif_each_gene_set get_mutation_count_for_peakmotif_nonpeakmotif_other_each_gene_set get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region summary_gene_ATCG_content get_merged_motif_bases_and_motif_count
# VARIABLE NAMING TODO:
export working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12_exclude_pseudo"
export motif_pattern="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/motifs/GAWGAW"
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export mutation_density_global_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/gene_6_sets
export variant_dir=$global_output_dir/$popu_symbol/variant_two_group

export genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
export gene_sets="DR HDR RNA_translation STEVOR VAR RIF"
export inconsistent_region_in_core_and_noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore/inconsistent_region_in_core_and_noncore.bed"
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export peak_loc=$m6mA_bam_dir/peak.bed
export seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
export parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
$parallel -j 6 get_mutation_count_for_motif_nonmotif_each_gene_set ::: $gene_sets ::: snps
$parallel -j 6 get_mutation_count_for_peakmotif_nonpeakmotif_other_each_gene_set ::: $gene_sets ::: snps
$parallel -j 6 summary_gene_ATCG_content ::: $gene_sets
$parallel -j 6 get_mutation_count_for_Mergedmotif_nonmotif_each_gene_set ::: $gene_sets ::: snps
$parallel -j 6 get_merged_motif_bases_and_motif_count ::: $gene_sets

##探究为什么信号区换成motif，就不好使了(TдT)
$parallel -j 6 explore_why ::: $gene_sets
