#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
make_window() {
    chrom=$1

    cd $window_bed_dir
    $bedtools makewindows -b $core_two_outgroup_dir/${chrom}.peak.bed -w $window_size -s $step_size -i srcwinnum >${chrom}.window.peak.bed
    $bedtools makewindows -b $core_two_outgroup_dir/${chrom}.nopeak.bed -w $window_size -s $step_size -i srcwinnum >${chrom}.window.nopeak.bed
}

split_window() {
    chrom=$1

    cd $window_bed_dir

    split -l 10000 ${chrom}.window.peak.bed -a 3 -d ${chrom}.window.peak.bed-
    split -l 10000 ${chrom}.window.nopeak.bed -a 3 -d ${chrom}.window.nopeak.bed-
}

transvert_motif_loc() {
    cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_GAWGAW/motif_loc
    for chrom in $(echo $autosomes | tr " " "\n"); do
        awk '{print $0"\t""motif_"NR}' ${chrom}.motif.bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$7}' >${chrom}.motif.right.bed
    done

    cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/motif_pattern_GAWGAW_nopeak/motif_loc
    for chrom in $(echo $autosomes | tr " " "\n"); do
        awk '{print $0"\t""motif_"NR}' ${chrom}.nopeak.motif.bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$7}' >${chrom}.nopeak.motif.right.bed
    done
}

submit_summary() {
    chrom=$1

    cd $window_bed_dir

    for prefix in $(echo "peak nopeak" | tr " " "\n"); do
        part_count=$(ls ${chrom}.window.${prefix}.bed-* | wc -l)
        let file_max_order=part_count-1
        for order in $(seq 0 $file_max_order); do
            sed -e "s/window_size_moban/$window_size/g" -e "s/step_size_moban/$step_size/g" -e "s/chrom_moban/$chrom/g" -e "s@output_dir_moban@$window_bed_dir@g" -e "s/order_moban/$order/g" -e "s/prefix_moban/$prefix/g" $code_dir/ss4_summary_for_each_part_in_chrom.sh >$code_dir/ss4_summary_for_each_part_in_chrom.sh_${window_size}_${step_size}_${chrom}_${prefix}_${order}
            sbatch $code_dir/ss4_summary_for_each_part_in_chrom.sh_${window_size}_${step_size}_${chrom}_${prefix}_${order}
        done
    done
}

integrate_and_clean() {
    cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_50_core_genome_two_outgroup
    rm -rf *window.*.bed-*
    tar -zcvf allChroms.window.original.summary.tar.gz *window.*.summary*

    for chrom in $(echo $autosomes | tr " " "\n"); do

        cat ${chrom}.window.peak.summary.* | sort -k1,1 -k2,2n >${chrom}.window.peak.summary
        cat ${chrom}.window.nopeak.summary.* | sort -k1,1 -k2,2n >${chrom}.window.nopeak.summary

        rm -rf ${chrom}.window.*.summary.*
    done

    ##这里的allChrom.window.peak.filter.summary跟allChroms.window.original.summary完全一样，程序里没有过滤
    cat *.window.peak.summary | sed -i 's/\t\+/\t/g' >allChrom.window.peak.filter.summary
    cat *.window.nopeak.summary | sed -i 's/\t\+/\t/g' >allChrom.window.nopeak.filter.summary

    sed -i 's/ /\t/g' allChrom.window.peak.filter.summary allChrom.window.nopeak.filter.summary
    gzip *summary *bed
}
export -f make_window split_window
# VARIABLE NAMING TODO:
export code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7/motif_Predict_peak"
export core_two_outgroup_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent"
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

export window_bed_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_50_core_genome_two_outgroup
export window_size=50
export step_size=10

export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
mkdir -p $window_bed_dir
transvert_motif_loc

$parallel make_window ::: $autosomes
$parallel split_window ::: $autosomes
$parallel submit_summary ::: $autosomes
integrate_and_clean
