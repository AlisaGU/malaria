#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err

# FUNCTIONS TODO:
get_submit_flank_in_two_outgroup_consistent_region() {
    mkdir -p $submit_flank_in_two_outgroup_consistent_dir
    for chrom in $(echo $autosomes | tr " " "\n"); do
        for por in $(echo $relative_por | tr " " "\n"); do
            mkdir -p $submit_flank_in_two_outgroup_consistent_dir/relative.${por}
            $bedtools intersect -a $submit_flank_whole_region_dir/relative.${por}/${chrom}.peak.relative.${por}.bed -b $two_outgroup_consistent_dir/${chrom}.peak.bed >$submit_flank_in_two_outgroup_consistent_dir/relative.${por}/${chrom}.peak.relative.${por}.bed
        done
        for len in $(echo $absolute_len | tr " " "\n"); do
            mkdir -p $submit_flank_in_two_outgroup_consistent_dir/absolute.${len}
            $bedtools intersect -a $submit_flank_whole_region_dir/absolute.${len}/${chrom}.peak.absolute.${len}.bed -b $two_outgroup_consistent_dir/${chrom}.peak.bed >$submit_flank_in_two_outgroup_consistent_dir/absolute.${len}/${chrom}.peak.absolute.${len}.bed
        done
    done
}

get_average_variant_count() {
    var_type=$1

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""

        regions="$two_outgroup_consistent_dir/${chrom}.nopeak.bed"
        for por in $(echo $relative_por | tr " " "\n"); do
            regions+=" $submit_flank_in_two_outgroup_consistent_dir/relative.${por}/${chrom}.peak.relative.${por}.bed"
        done
        for len in $(echo $absolute_len | tr " " "\n"); do
            regions+=" $submit_flank_in_two_outgroup_consistent_dir/absolute.${len}/${chrom}.peak.absolute.${len}.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            if [ $var_type == "transition" ] || [ $var_type == "transversion" ] || [ $var_type == "allvar" ] || [ $var_type == "snps" ] || [ $var_type == "indels" ]; then
                chrom_var_count+=" $(get_average_variant_count_for_chrom_whole_region $region_bed $popu_symbol)"
            else
                chrom_var_count+=" $(get_average_variant_count_for_chrom_base_count $region_bed $popu_symbol $var_type)"
            fi
        done
        echo $chrom_var_count >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom_whole_region() {
    region_bed=$1
    popu_symbol=$2
    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}

get_average_variant_count_for_chrom_base_count() {
    region_bed=$1
    popu_symbol=$2
    var_type=$3

    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    ori_base=$(echo $var_type | awk -F"2" '{print $1}')
    chrom_region_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${ori_base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

    chrom_mutation_count=""
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6+$9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $6}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $9}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $5}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "insertion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $7}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "deletion" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $8}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $10}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $11}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $12}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $13}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $14}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "A2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $15}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2C" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $16}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "G2T" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $17}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $18}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "C2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $19}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2A" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $20}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "T2G" ]; then
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '{print $21}' | awk '{s+=$1} END {print s}')
    elif [ $var_type == "fixation" ]; then # 这里指snp的固定
        chrom_mutation_count=$($bedtools intersect -a $variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed -b $region_bed | awk '$6==1 && $9==0{print $0}' | wc -l)
    fi

    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi

    chrom_region_mutation_count_mean=$(echo "scale=10;${chrom_mutation_count}/${chrom_region_len}" | bc)
    chrom_region_mutation_count_mean1=$(printf "%0.6f" $chrom_region_mutation_count_mean)
    echo $chrom_region_mutation_count_mean1
}

compute_base_proporation() {
    base=$1
    base_proporation=$base_proporation_dir/${base}_mean_proporation
    rm -rf $base_proporation

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_base_propor=""

        regions="$two_outgroup_consistent_dir/${chrom}.nopeak.bed"
        for por in $(echo $relative_por | tr " " "\n"); do
            regions+=" $submit_flank_in_two_outgroup_consistent_dir/relative.${por}/${chrom}.peak.relative.${por}.bed"
        done
        for len in $(echo $absolute_len | tr " " "\n"); do
            regions+=" $submit_flank_in_two_outgroup_consistent_dir/absolute.${len}/${chrom}.peak.absolute.${len}.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            whole_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
            base_len=$($seqkit subseq --quiet --bed $region_bed $genome | $seqkit fx2tab -i -H -n -C "${base}" | grep -v "^#" | awk '{print $2}' | awk '{s+=$1} END {print s}')

            base_proportion=$(echo "scale=10;${base_len}/${whole_region_len}" | bc)
            base_proportion1=$(printf "%0.6f" $base_proportion)
            chrom_base_propor+=" $(echo $base_proportion1)"
        done

        echo $chrom_base_propor >>$base_proporation
    done
}
# VARIABLE NAMING TODO:
popu_symbol="WSEA"
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output_jiang"
submit_flank_whole_region_dir=${macs2_out_dir}/chrom
two_outgroup_consistent_dir=${macs2_out_dir}/two_outgroup_consistent
submit_flank_in_two_outgroup_consistent_dir=$macs2_out_dir/submit_flank_in_3D7_african_strain_consistent_dir
variant_dir=$global_output_dir/$popu_symbol/variant_two_group
mutation_density_dir=$global_output_dir/$popu_symbol/mutation_density_two_group_consistent_two_group_consistent
base_proporation_dir=$global_output_dir/$popu_symbol/base_proporation_two_group_consistent_two_group_consistent

absolute_len="10 30 50 100 200"
relative_por="0.01 0.05 0.1 0.2 0.5"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
mkdir -p $mutation_density_dir $base_proporation_dir

# get_submit_flank_in_two_outgroup_consistent_region # OCE群体已经做过这一步，不需要重复做

get_average_variant_count "allvar"
get_average_variant_count "snps"
get_average_variant_count "indels"
# get_average_variant_count "transition"
# get_average_variant_count "transversion"
# $code_dir/s2_get_ts_tv_ratio_each_popu.R $popu_symbol
get_average_variant_count "A2G"
get_average_variant_count "G2A"
get_average_variant_count "T2C"
get_average_variant_count "C2T"
get_average_variant_count "A2C"
get_average_variant_count "A2T"
get_average_variant_count "G2C"
get_average_variant_count "G2T"
get_average_variant_count "C2A"
get_average_variant_count "C2G"
get_average_variant_count "T2A"
get_average_variant_count "T2G"

compute_base_proporation "A"
compute_base_proporation "T"
compute_base_proporation "C"
compute_base_proporation "G"
