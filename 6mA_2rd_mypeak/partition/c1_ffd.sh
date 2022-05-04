#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
check_cds_pep_seqid_identical() {
    protein_md5sum_id=$(grep "^>" ${pep_seq} | awk '{print $1}' | sed 's/>//g' | sed 's/-p1//g' | md5sum | awk '{print $1}')
    cds_md5sum_id=$(grep "^>" ${cds_seq} | awk '{print $1}' | sed 's/>//g' | md5sum | awk '{print $1}')
    if [ $protein_md5sum_id == $cds_md5sum_id ]; then
        echo "sequence id of protein and cds is identical"
    else
        echo "NO"
    fi
}

get_ffd_bed_for_one_mRNA() {
    id=$1
    # 检查cds直接翻译是否和网站提供的pep序列一样
    cds_translate_status=$(check_cds_translate $id)

    if [ $cds_translate_status -eq "1" ]; then
        mRNA_cds_seq=$output_dir/$id".cds.fa"
        $seqkit grep -p $id $cds_seq >$mRNA_cds_seq
        $code_dir/s1_get_ffd_bed.R $id $output_dir
        rm -rf $mRNA_cds_seq
    fi
}

check_cds_translate() {
    id=$1
    cds_translate_md5sum=$(seqkit grep -p $id $cds_seq | seqkit translate | grep -v "^>" | sed '$ s/*$//' | grep -v "^$" | md5sum | awk '{print $1}')
    pep_md5sum=$(seqkit grep -p $id"-p1" $pep_seq | grep -v "^>" | md5sum | awk '{print $1}')
    if [ $cds_translate_md5sum == $pep_md5sum ]; then
        echo "1"
    else
        echo "0"
    fi
}

get_all_ffd_site_bed() {
    grep "^>" ${cds_seq} | awk '{print $1}' | sed 's/>//g' >$all_mRNA_id

    for i in $(cat $all_mRNA_id); do
        get_ffd_bed_for_one_mRNA $i
    done

    # 注意科学计数法
    cat $output_dir/*bed | sort -k1,1 -k2,2n >$ffd_bed
    rm -rf $output_dir/*ffd.bed
}

select_peak_nopeak_ffd_site_bed() {
    mkdir -p $ffd_bed_dir
    for autosome in $(echo $autosomes | tr " " "\n"); do
        $bedtools intersect -a $ffd_bed -b $peak_nopeak_bed_dir/$autosome.nopeak.bed >$ffd_bed_dir/${autosome}.nopeak.ffd.bed
    done

    for absolute in $(echo $absolute_len | tr " " "\n"); do
        specific_peak_flank_dir=$ffd_bed_dir/"absolute."$absolute
        mkdir -p $specific_peak_flank_dir
        for autosome in $(echo $autosomes | tr " " "\n"); do
            $bedtools intersect -a $ffd_bed -b $peak_nopeak_bed_dir/absolute.$absolute/$autosome.* >$specific_peak_flank_dir/${autosome}.peak.ffd.bed
        done
    done

    for relative in $(echo $relative_por | tr " " "\n"); do
        specific_peak_flank_dir=$ffd_bed_dir/"relative."$relative
        mkdir -p $specific_peak_flank_dir
        for autosome in $(echo $autosomes | tr " " "\n"); do
            $bedtools intersect -a $ffd_bed -b $peak_nopeak_bed_dir/relative.$relative/$autosome.* >$specific_peak_flank_dir/${autosome}.peak.ffd.bed
        done
    done
}

get_average_variant_count() {
    for popu_symbol in $(echo $popu_symbols); do
        get_average_variant_count_each_popu $popu_symbol "allvar"
        get_average_variant_count_each_popu $popu_symbol "snps"
        get_average_variant_count_each_popu $popu_symbol "indels"
        # get_average_variant_count_each_popu $popu_symbol "transition"
        # get_average_variant_count_each_popu $popu_symbol "transversion"
        # $code_dir/s1_get_ts_tv_ratio_each_popu.R $popu_symbol
    done
}

get_average_variant_count_each_popu() {
    popu_symbol=$1
    var_type=$2

    mutation_density_dir=$output_dir/$popu_symbol/mutation_density
    mkdir -p $mutation_density_dir

    mean_variant_count=$mutation_density_dir/${var_type}_mean_variant_count
    rm -rf $mean_variant_count

    for chrom in $(echo $autosomes | tr " " "\n"); do
        chrom_var_count=""

        regions="$ffd_bed_dir/${chrom}.nopeak.ffd.bed"
        for por in $(echo $relative_por | tr " " "\n"); do
            regions+=" $ffd_bed_dir/relative.${por}/${chrom}.peak.ffd.bed"
        done
        for len in $(echo $absolute_len | tr " " "\n"); do
            regions+=" $ffd_bed_dir/absolute.${len}/${chrom}.peak.ffd.bed"
        done

        for region_bed in $(echo $regions | tr " " "\n"); do
            chrom_var_count+=" $(get_average_variant_count_for_chrom $region_bed $popu_symbol)"
        done
        echo $chrom_var_count >>$mean_variant_count
    done
    # $code_dir/s6_plot.R $mean_variant_count
}

get_average_variant_count_for_chrom() {
    region_bed=$1
    popu_symbol=$2
    chrom_number=$(basename $region_bed | grep -o "Pf3D7_.*_v3" | sed 's/Pf3D7_//g' | sed 's/_v3//g')

    chrom_mutation_count_bed=""
    variant_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/1_2rd_initial_evaluation/${popu_symbol}/variant"
    if [ $var_type == "allvar" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_polysite.bed
    elif [ $var_type == "snps" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_snp_polysite.bed
    elif [ $var_type == "indels" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_indel_polysite.bed
    elif [ $var_type == "transition" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_ts_polysite.bed
    elif [ $var_type == "transversion" ]; then
        chrom_mutation_count_bed=$variant_dir/${popu_symbol}_PASS_${chrom_number}_tv_polysite.bed
    fi

    chrom_region_len=$(awk '{print $3-$2}' $region_bed | awk '{s+=$1} END {print s}')
    chrom_mutation_count=$($bedtools intersect -a $chrom_mutation_count_bed -b $region_bed | awk '{print $4}' | awk '{s+=$1} END {print s}')
    if [ $chrom_mutation_count -gt 0 ] 2>/dev/null; then
        echo ""
    else
        chrom_mutation_count=0
    fi
    chrom_region_mutation_count_mean=$(echo "scale=4;${chrom_mutation_count}/${chrom_region_len}" | bc)
    echo $chrom_region_mutation_count_mean
}

plot() {
    $code_dir/s1_plot_6_population_ts_tv_mutation_density.R
}
# VARIABLE NAMING TODO:
code_dir="/picb/evolgen/users/gushanshan/projects/malaria/code/6mA_2rd/partition"
macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/2_2rd_partition/ffd"
global_seq_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome"
peak_nopeak_bed_dir=$macs2_out_dir/chrom
ffd_bed_dir=$macs2_out_dir/ffd

cds_seq=$global_seq_dir/"PlasmoDB-36_Pfalciparum3D7_AnnotatedCDSs.fasta"
pep_seq=$global_seq_dir/"PlasmoDB-36_Pfalciparum3D7_AnnotatedProteins.fasta"
all_mRNA_id=$output_dir/"all_mRNA_id"
absolute_len="10 30 50 100 200"
relative_por="0.01 0.05 0.1 0.2 0.5"
ffd_bed=$output_dir/"four_fold_degenerate_site.bed"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
popu_symbols="CAF EAF WAF WAF_CAF_EAF WSEA SAM_WAF_CAF_EAF_SAS_WSEA_ESEA_OCE"

seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -x
# check_cds_pep_seqid_identical
# get_all_ffd_site_bed
# select_peak_nopeak_ffd_site_bed
get_average_variant_count
# plot
