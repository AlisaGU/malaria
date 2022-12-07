#!/bin/bash
#SBATCH -n 60
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_sliding_window_for_chrom() {
    chrom=$1
    flank_length_type="relative"
    flank_length="0.25"
    peak_pos=$chrom_bed_dir/${chrom}.peak.bed
    treatment_level=$modifi_level_dir/${chrom}.treat.bed
    control_level=$modifi_level_dir/${chrom}.cont.bed

    output_dir=$chrom_bed_dir/sliding_window/${flank_length_type}.${flank_length}/${chrom}
    mkdir -p $output_dir
    $code_dir/s0_get_submit_flank25Perc.R $chrom $peak_pos $flank_length_type $flank_length $treatment_level $control_level $output_dir
}

export -f get_sliding_window_for_chrom
# VARIABLE NAMING TODO:
export code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7/gene_groups_analysis/whole_genome_peak25Perc"
export macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
export chrom_bed_dir=$macs2_out_dir/chrom_wholeGenome
export modifi_level_dir=$macs2_out_dir/modifi_level ## modifi_level文件夹创建的时候，就没有核心不核心的概念

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
$parallel -j 14 get_sliding_window_for_chrom ::: $autosomes

cd $chrom_bed_dir/sliding_window/relative.0.25
cat $(find ./ -name *.bed) | sort -k1,1 -k2,2n >sliding_window_25Perc.bed
