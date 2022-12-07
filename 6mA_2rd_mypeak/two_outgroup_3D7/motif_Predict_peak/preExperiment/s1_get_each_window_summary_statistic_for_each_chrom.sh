#!/bin/bash
#SBATCH -n 1
#SBATCH -p mqueue
#SBATCH -w mnode019
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
make_window() {
    grep $chrom $genome_length | $bedtools makewindows -g stdin -w $window_size -s $step_size -i srcwinnum >${chrom}.window.bed
}

split_bdg() {
    mkdir -p $bdg_dir
    cd $bdg_dir
    sort -k1,1 -k2,2n $treat_bdg >${treat_bdg}.sort
    sort -k1,1 -k2,2n $cont_bdg >${cont_bdg}.sort
    awk '{print >$1".treat.bdg"}' ${treat_bdg}.sort
    awk '{print >$1".cont.bdg"}' ${cont_bdg}.sort
    rm -rf ${treat_bdg}.sort ${cont_bdg}.sort
    cd -
}

split_window() {
    split -l 10000 ${chrom}.window.bed -a 3 -d ${chrom}.window.bed-
}

submit_summary() {
    part_count=$(ls ${chrom}.window.bed-* | wc -l)
    let file_max_order=part_count-1
    for order in $(seq 0 $file_max_order); do
        sed -e "s/window_size_moban/$window_size/g" -e "s/step_size_moban/$step_size/g" -e "s/chrom_moban/$chrom/g" -e "s@output_dir_moban@$output_dir@g" -e "s/order_moban/$order/g" $code_dir/ss1_summary_for_each_part_in_chrom.sh >$code_dir/ss1_summary_for_each_part_in_chrom.sh_${window_size}_${step_size}_${chrom}_${order}
        sbatch $code_dir/ss1_summary_for_each_part_in_chrom.sh_${window_size}_${step_size}_${chrom}_${order}
    done
}
# VARIABLE NAMING TODO:
chrom=chrom_template
window_size=window_size_template ## 单位：bp
step_size=step_size_template     ## 单位：bp
output_dir=output_dir_template

code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7/motif_Predict_peak"

motif_genome_loc=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/motif_GAWGAW_genome_loc
genome_length=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/chrom.length
treat_bdg=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7-T3_treat_pileup.bdg
cont_bdg=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/3D7-T3_control_lambda.bdg
whole_genome_peak_bed=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed
bdg_dir=$output_dir/../../whole_genome_bdg
seqkit=/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
datamash=/picb/evolgen/users/gushanshan/software/datamash/bin/datamash
# VARIABLE NAMING for test module TODO:
# chrom="Pf3D7_01_v3"
# window_size=100
# step_size=10
# output_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/${window_size}_${step_size}
# PROCESS TODO:
# set -ex
# split_bdg

mkdir -p $output_dir
cd $output_dir
# make_window
# split_window

##统计各个窗口
submit_summary
