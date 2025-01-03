#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
locate_motif_in_peak_and_control_two_outgroup_consistent() {
    mkdir -p $motif_bed_dir
    cd $motif_bed_dir

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $core_peak_seq_dir/${autosome}.peak.core.consistent.fa | $seqkit locate -i -d -f $motif --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6"\t"$4}' | sort -k1,1 -k2,2n >${autosome}.motif.peak.bed
    done

    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $core_nopeak_seq_dir/${autosome}.nopeak.core.consistent.fa | $seqkit locate -i -d -f $motif --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6"\t"$4}' | sort -k1,1 -k2,2n >${autosome}.motif.nopeak.bed
    done
    cd -
}

# VARIABLE NAMING TODO:
# motif=motif_template
code_dir="/picb/evolgen2/users/gushanshan/projects/malaria/code/6mA_2rd/two_outgroup_3D7"
motif_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/motifs"
motif_seq=motif_seq_template
motif=$motif_dir/$motif_seq

macs2_out_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output"
macs2_two_outgroup_consistent_dir=$macs2_out_dir/two_outgroup_consistent
single_motif_global_output_dir=$macs2_two_outgroup_consistent_dir/single_motif_pattern
single_motif_output_dir=$single_motif_global_output_dir/$motif_seq
motif_bed_dir=$single_motif_output_dir/motif_bed

core_peak_seq_dir=$macs2_two_outgroup_consistent_dir/peak_seq
core_nopeak_seq_dir=$macs2_two_outgroup_consistent_dir/nopeak_seq

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
mkdir -p $single_motif_output_dir
locate_motif_in_peak_and_control_two_outgroup_consistent
