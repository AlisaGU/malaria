#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_whole_genome_1peak_2nopeak_3nopeak_away_from_peak() {
    mkdir -p $chrom_bed_wholeGenome_dir
    cd $chrom_bed_wholeGenome_dir
    awk '{print >$1".peak.bed"}' $peak
    awk '{print >$1".nopeak.bed"}' $nopeak
    awk '{print $0"\t"$3-$2}' $nopeak | awk '$4>1000{s=$2+500;e=$3-500;print $1"\t"s"\t"e >$1".nopeak.awayFromPeak.bed"}'
}

make_window() {
    window_size=$1
    chrom=$2

    window_bed_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_${window_size}
    mkdir -p $window_bed_dir
    cd $window_bed_dir
    $bedtools makewindows -b $chrom_bed_wholeGenome_dir/${chrom}.peak.bed -w $window_size -s $step_size -i srcwinnum >${chrom}.window.peak.bed
    $bedtools makewindows -b $chrom_bed_wholeGenome_dir/${chrom}.nopeak.awayFromPeak.bed -w $window_size -s $step_size -i srcwinnum >${chrom}.window.nopeak.awayFromPeak.bed
}

split_window() {
    window_size=$1
    chrom=$2
    window_bed_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_${window_size}
    mkdir -p $window_bed_dir
    cd $window_bed_dir

    split -l 10000 ${chrom}.window.peak.bed -a 3 -d ${chrom}.window.peak.bed-
    split -l 10000 ${chrom}.window.nopeak.awayFromPeak.bed -a 3 -d ${chrom}.window.nopeak.awayFromPeak.bed-
}

submit_summary() {
    window_size=$1
    chrom=$2
    window_bed_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_${window_size}
    mkdir -p $window_bed_dir
    cd $window_bed_dir

    for prefix in $(echo "peak nopeak.awayFromPeak" | tr " " "\n"); do
        part_count=$(ls ${chrom}.window.${prefix}.bed-* | wc -l)
        let file_max_order=part_count-1
        for order in $(seq 0 $file_max_order); do
            sed -e "s/window_size_moban/$window_size/g" -e "s/step_size_moban/$step_size/g" -e "s/chrom_moban/$chrom/g" -e "s@output_dir_moban@$window_bed_dir@g" -e "s/order_moban/$order/g" -e "s/prefix_moban/$prefix/g" $code_dir/ss1_summary_for_each_part_in_chrom.sh >$code_dir/ss1_summary_for_each_part_in_chrom.sh_${window_size}_${step_size}_${chrom}_${prefix}_${order}
            sbatch $code_dir/ss1_summary_for_each_part_in_chrom.sh_${window_size}_${step_size}_${chrom}_${prefix}_${order}
        done
    done
}

integrate_and_clean() {
    window_size=$1
    cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/peak_nopeak_${window_size}
    rm -rf *window.*.bed-*
    tar -zcvf allChroms.window.original.summary.tar.gz *window.*.summary*

    for chrom in $(echo $autosomes | tr " " "\n"); do

        cat ${chrom}.window.peak.summary.* | awk 'NF==8' | sort -k1,1 -k2,2n >${chrom}.window.peak.summary
        cat ${chrom}.window.nopeak.awayFromPeak.summary.* | awk 'NF==8' | sort -k1,1 -k2,2n >${chrom}.window.nopeak.awayFromPeak.summary

        rm -rf ${chrom}.window.*.summary.*
    done

    cat *.window.peak.summary >allChrom.window.peak.filter.summary                               ##filter指删掉了没有fold enrichment的窗口。这种窗口一般位于染色体末端，没有对应的bedgraph数据
    cat *.window.nopeak.awayFromPeak.summary >allChrom.window.nopeak.awayFromPeak.filter.summary ##filter指删掉了没有fold enrichment的窗口。这种窗口一般位于染色体末端，没有对应的bedgraph数据

    sed -i 's/ /\t/g' allChrom.window.peak.filter.summary allChrom.window.nopeak.awayFromPeak.filter.summary
    gzip *summary *bed
}
export -f make_window split_window submit_summary
# VARIABLE NAMING TODO:
# window_sizes="100 200 300 400 500"
window_sizes="50"

export step_size=10

export code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7/motif_Predict_peak"
macs2_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
peak=$macs2_output_dir/peak.bed
nopeak=$macs2_output_dir/nopeak.bed

export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

export chrom_bed_wholeGenome_dir=$macs2_output_dir/chrom_wholeGenome
export bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
get_whole_genome_1peak_2nopeak_3nopeak_away_from_peak
$parallel make_window ::: $window_sizes ::: $autosomes
$parallel split_window ::: $window_sizes ::: $autosomes
$parallel submit_summary ::: "50" ::: $autosomes
integrate_and_clean
