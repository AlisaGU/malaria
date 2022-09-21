#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_gene_locus() {
    gene_set=$1
    cd $gene_sets_global_dir/$gene_set/
    grep -f ${gene_set}_genes_list $anno | awk '$3=="gene"' | awk '{split($9,a,";");split(a[1],b,"=");print $1"\t"$4-1"\t"$5"\t"b[2]}' | sort -k1,1 -k2,2n >${gene_set}_genes.bed
}

get_mutation_count_for_each_gene_set() {
    gene_set=$1
    var_type=$2

    mutation_density_dir=$mutation_density_global_dir/${gene_set}
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    cat $gene_sets_global_dir/$gene_set/${gene_set}_genes.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $peak_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA" "NA" "NA")"
        elif [ $peak_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            #peak region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+="$(echo $gene_name "$result")"
            #control region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_loc | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
        fi
        echo $value >>$mean_variant_count
    done
}

get_mutation_count_for_each_gene_set_include_nopeak_gene() {
    gene_set=$1
    var_type=$2

    mutation_density_dir=$mutation_density_global_dir/${gene_set}
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    cat $gene_sets_global_dir/$gene_set/${gene_set}_genes.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $peak_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA")"

            chrom=$(echo $line | awk '{print $1}')

            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_loc | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
        elif [ $peak_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            #peak region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_loc | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+="$(echo $gene_name "$result")"
            #control region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_loc | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result")"
        fi
        echo $value >>$mean_variant_count
    done
}

get_mutation_count_for_each_gene_set_include_nopeak_gene_with_fold_enrichment() {
    gene_set=$1
    var_type=$2

    mutation_density_dir=$mutation_density_global_dir/${gene_set}
    mkdir -p $mutation_density_dir
    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    cat $gene_sets_global_dir/$gene_set/${gene_set}_genes.bed | while read line; do
        peak_exist_in_gene_status=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | wc -l)
        value=""
        gene_name=$(echo $line | awk '{print $4}')
        if [ $peak_exist_in_gene_status -eq 0 ]; then
            value+="$(echo $gene_name "NA" "NA")"

            chrom=$(echo $line | awk '{print $1}')

            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            value+=" $(echo "$result" "NA")"
        elif [ $peak_exist_in_gene_status -gt 0 ]; then
            chrom=$(echo $line | awk '{print $1}')
            #peak region
            bed=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")
            mean_fold_enrichment=$(echo $"$line" | $bedtools intersect -a stdin -b $peak_with_fold_enrichment -wb | awk '{print $NF}' | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }')
            value+="$(echo $gene_name "$result")"
            #control region
            bed=$(echo $"$line" | $bedtools subtract -a stdin -b $peak_with_fold_enrichment | $bedtools subtract -a stdin -b $inconsistent_region_in_core_and_noncore)
            result=$(get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region $chrom $var_type "$bed")

            value+=" $(echo "$result" $mean_fold_enrichment)"
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
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $22}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $23}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $24}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $25}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $26}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $27}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $28}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_single_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $29}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_insert" ]; then # 这里指snp的固定
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $30}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $31}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $32}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_insert" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $33}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_multi_del" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $37}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $22+$26+$30+$34}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $23+$27+$31+$35}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $24+$28+$32+$36}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G_indel" ]; then
        chrom_mutation_count=$(echo $"$region_bed" | $bedtools intersect -a $variant_dir/${popu_symbol}_VQSLOD_gt0_${chrom_number}_polysite.bed -b stdin | awk '{print $25+$29+$33+$37}' | awk '{s+=$1} END {print s}')
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        # In this case, the ":" is just a "do-nothing" place holder
        :
    else
        chrom_mutation_count=0
    fi

    echo "$chrom_mutation_count $chrom_region_len"
}
export -f get_gene_locus get_mutation_count_for_each_gene_set get_average_variant_count_for_chrom_two_outgroup_consistent_in_core_and_noncore_whole_region get_mutation_count_for_each_gene_set_include_nopeak_gene get_mutation_count_for_each_gene_set_include_nopeak_gene_with_fold_enrichment
# VARIABLE NAMING TODO:
export popu_symbol="ESEA_WSEA_OCE_SAM_SAS"
export global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
export m6mA_bam_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export peak_loc=$m6mA_bam_dir/peak.bed
export peak_with_fold_enrichment=$m6mA_bam_dir/peak_with_fold_enrichment.bed
export mutation_density_global_dir=$global_output_dir/${popu_symbol}/two_outgroup_consistent_in_core_and_noncore/gene_5_sets
export variant_dir=$global_output_dir/$popu_symbol/variant_two_group

export gene_sets_global_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/gene_groups/gene_5_sets_2022_07_12"
export inconsistent_region_in_core_and_noncore="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/3D7_distant_africa_strain_consistent_region_in_core_and_noncore/inconsistent_region_in_core_and_noncore.bed"
export anno="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/PlasmoDB-36_Pfalciparum3D7.gff"
gene_sets="DR HDR RNA_translation STEVOR VAR RIF" ## 跟c11_gene_5_sets相比，只是加了一个RIF家族
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
$parallel get_gene_locus ::: $gene_sets
$parallel -j 5 get_mutation_count_for_each_gene_set ::: $gene_sets ::: snps
$parallel -j 5 get_mutation_count_for_each_gene_set_include_nopeak_gene ::: $gene_sets ::: snps
$parallel -j 5 get_mutation_count_for_each_gene_set_include_nopeak_gene_with_fold_enrichment ::: $gene_sets ::: snps
