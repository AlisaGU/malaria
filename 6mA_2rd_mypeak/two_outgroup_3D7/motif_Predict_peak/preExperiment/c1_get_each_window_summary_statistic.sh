#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
locate_GAWGAW_in_genome() {
    cd $motif_genome_loc
    cat $genome | $seqkit locate -i -d -f $motif_pattern --bed | awk '{print >$1".genome.motif.bed"}'
}

summary_each_chrom_window_statistic() {
    for window_size in $(echo $window_sizes | tr " " "\n"); do
        for step_size in $(echo $step_sizes | tr " " "\n"); do
            for chrom in $(echo $autosomes | tr " " "\n" | awk 'NR>=9'); do
                output_dir=$result_dir/${window_size}_${step_size}
                sed -e "s/window_size_template/$window_size/g" -e "s/step_size_template/$step_size/g" -e "s/chrom_template/$chrom/g" -e "s@output_dir_template@$output_dir@g" $code_dir/s1_get_each_window_summary_statistic_for_each_chrom.sh >$code_dir/s1_get_each_window_summary_statistic_for_each_chrom.sh_${window_size}_${step_size}_${chrom}
                sbatch $code_dir/s1_get_each_window_summary_statistic_for_each_chrom.sh_${window_size}_${step_size}_${chrom}
            done
        done
    done
}

integrate_and_clean() {
    cd $result_dir/${window_size}_${step_size}
    rm -rf *window.bed-*
    tar -zcvf allChroms.window.original.summary.tar.gz *window.summary*

    for chrom in $(echo $autosomes | tr " " "\n"); do
        cat ${chrom}.window.summary.* | awk 'NF==9' | sort -k1,1 -k2,2n >${chrom}.window.summary
        rm -rf ${chrom}.window.summary.*
    done

    cat *.window.summary >allChrom.window.filter.summary ##filter指删掉了没有fold enrichment的窗口。这种窗口一般位于染色体末端，没有对应的bedgraph数据
    sed -i 's/ /\t/g' allChrom.window.filter.summary
    gzip *summary *bed
}
# VARIABLE NAMING TODO:
window_sizes=300 ## 单位：bp
step_sizes=20    ## 单位：bp

code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd_mypeak/two_outgroup_3D7/motif_Predict_peak"
result_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW
peak_nopeak_bed_dir=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom
motif_genome_loc=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/motif_GAWGAW_genome_loc

export motif_pattern="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/motifs/GAWGAW"
genome=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

seqkit=/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
# locate_GAWGAW_in_genome
summary_each_chrom_window_statistic
# integrate_and_clean
