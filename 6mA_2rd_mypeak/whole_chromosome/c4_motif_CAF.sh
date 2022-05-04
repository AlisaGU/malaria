#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
extract_peak_seq_in_core_region() {
    mkdir -p $core_peak_seq_dir
    cd $core_peak_seq_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        $seqkit subseq --bed $core_peak_nopeak_chrom_dir/${autosome}.peak.bed $genome -o ${autosome}.peak.core.fa
    done

    cd -
}

locate_motif() {
    mkdir -p $motif_bed_dir
    cd $motif_bed_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $core_peak_seq_dir/${autosome}.peak.core.fa | $seqkit locate -i -d -f $motif_pattern --bed >${autosome}.motif.bed
    done

    cd -
}

get_block_bed_for_motif() {
    mkdir -p $region_bed_dir
    for bed in $(ls $motif_bed_dir); do
        $code_dir/s4_get_block_bed_for_motif.R ${motif_bed_dir}/$bed $region_bed_dir
    done
}

get_average_variant_count() {
    for popu_symbol in $(echo $popu_symbols); do
        get_average_variant_count_each_popu $popu_symbol "transition"
        get_average_variant_count_each_popu $popu_symbol "transversion"
        # $code_dir/s2_get_ts_tv_ratio_each_popu.R $popu_symbol
        get_average_variant_count_each_popu $popu_symbol "A2G"
        get_average_variant_count_each_popu $popu_symbol "G2A"
        get_average_variant_count_each_popu $popu_symbol "T2C"
        get_average_variant_count_each_popu $popu_symbol "C2T"
        get_average_variant_count_each_popu $popu_symbol "A2C"
        get_average_variant_count_each_popu $popu_symbol "A2T"
        get_average_variant_count_each_popu $popu_symbol "G2C"
        get_average_variant_count_each_popu $popu_symbol "G2T"
        get_average_variant_count_each_popu $popu_symbol "C2A"
        get_average_variant_count_each_popu $popu_symbol "C2G"
        get_average_variant_count_each_popu $popu_symbol "T2A"
        get_average_variant_count_each_popu $popu_symbol "T2G"
        get_average_variant_count_each_popu $popu_symbol "fixation"
        # $code_dir/s2_compare_mutation.R $popu_symbol
    done
}

get_average_variant_count_each_popu() {
    popu_symbol=$1
    var_type=$2

    mutation_density_dir=$working_dir/mutation_density
    mkdir -p $mutation_density_dir

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    # mean_variant_count=$mutation_density_dir/${var_type}_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""

        regions="$core_peak_nopeak_chrom_dir/${chrom}.nopeak.bed"
        for por in $(echo $regions_type | tr " " "\n"); do
            regions+=" $region_bed_dir/${chrom}_${por}.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom $region_bed $popu_symbol)"
        done
        echo $chrom_var_count >>$mean_variant_count
    done
}

get_average_variant_count_for_chrom() {
    region_bed=$1
    popu_symbol=$2
    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')
    variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}/variant"

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

    # echo $chrom_mutation_count
}
# VARIABLE NAMING TODO:
popu_symbols="CAF"
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/whole_chromosome"
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/"$popu_symbols"/motif"
# intergenic_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/intergenic/intergenic"
# mutation_density_global_dir=""
# core_peak_nopeak_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"

# core_peak_bed=$core_nopeak_bed_dir/"peak.core.bed"
core_peak_nopeak_chrom_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom"
core_peak_seq_dir=$working_dir/"peak_seq"
# core_peak_chrom_dir=
motif_bed_dir=$working_dir/motif_loc
motif_pattern="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition_jiang/intergenic_motif/motif_pattern"
region_bed_dir=$working_dir/region_loc

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
regions_type="collapsed_motif region1 region2 region3 region4"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
csvtk="/home/gushanshan/bin/csvtk"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# extract_peak_seq_in_core_region
# locate_motif
# get_block_bed_for_motif
get_average_variant_count
