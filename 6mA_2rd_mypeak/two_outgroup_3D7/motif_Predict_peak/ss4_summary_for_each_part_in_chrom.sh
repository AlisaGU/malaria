#!/bin/bash
#SBATCH -n 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
summary() {
    file=$(ls | grep ${chrom}.window.${prefix}.bed-[0]*$order$)
    rm -rf ${chrom}.window.${prefix}.summary.${order}

    motif_genome_loc=""
    if [ $prefix == "peak" ]; then
        motif_genome_loc=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_GAWGAW/motif_loc
    else
        motif_genome_loc=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_pattern_GAWGAW_nopeak/motif_loc
    fi

    cat $file | while read line; do
        ## 窗口信息
        let window_start=$(echo $line | awk '{print $2}')+1
        window_end=$(echo $line | awk '{print $3}')
        let window_start_minus_1=window_start-1

        ## 计算窗口里有多少碱基位于motif里
        base_count_in_motif=$(echo "$line" | $bedtools intersect -a stdin -b $motif_genome_loc/${chrom}.*motif.right.bed | sort -k1,1 -k2,2n | $bedtools merge -i stdin | awk '{s+=$3-$2}END{print s}')
        if [[ $base_count_in_motif -gt 0 ]]; then
            :
        else
            base_count_in_motif=0
        fi

        ## 计算窗口里突变count
        snp_mutation_count=$(get_average_variant_count_each_popu "$line" "snps" $prefix)
        indel_mutation_count=$(get_average_variant_count_each_popu "$line" "indels" $prefix)
        ##总结
        echo -e $chrom"\t"$window_start_minus_1"\t"$window_end"\t"$base_count_in_motif"\t"$snp_mutation_count"\t"$indel_mutation_count >>${chrom}.window.${prefix}.summary.${order}
    done

    rm -rf ${chrom}.${prefix}.${order}.*
}

get_average_variant_count_each_popu() {
    region_bed=$1
    var_type=$2
    prefix=$3
    chrom_number=$(echo $region_bed | awk '{print $1}' | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    if [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo "$region_bed" | $bedtools intersect -a $variant_dir/ESEA_WSEA_OCE_SAM_SAS_PASS_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo "$region_bed" | $bedtools intersect -a $variant_dir/ESEA_WSEA_OCE_SAM_SAS_PASS_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    echo $chrom_mutation_count
}
# VARIABLE NAMING TODO:
order=order_moban
chrom=chrom_moban
window_size=window_size_moban ## 单位：bp
step_size=step_size_moban     ## 单位：bp
output_dir=output_dir_moban
prefix=prefix_moban

variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/ESEA_WSEA_OCE_SAM_SAS/variant_two_group"

seqkit=/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
cd $output_dir
summary
