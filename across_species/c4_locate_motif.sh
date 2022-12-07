#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_peak_nopeak_in_two_outgroup_consistent_region() {
    mkdir -p $peak_nopeak_dir
    cd $peak_nopeak_dir
    for autosome in $(echo $autosomes | tr " " "\n"); do
        $bedtools intersect -a $two_outgroup_consistent_region_bed_dir/${autosome}.two_consistent_region.bed -b $m6mA_peak_bed >${autosome}.two_consistent_region.peak.bed
        $bedtools subtract -a $two_outgroup_consistent_region_bed_dir/${autosome}.two_consistent_region.bed -b $m6mA_peak_bed >${autosome}.two_consistent_region.nopeak.bed
    done
}
get_seq() {
    mkdir -p $seq_dir
    cd $seq_dir
    for autosome in $(echo $autosomes | tr " " "\n"); do
        $seqkit subseq --bed $peak_nopeak_dir/${autosome}.two_consistent_region.peak.bed $genome -o ${autosome}.two_consistent_region.peak.fa
        $seqkit subseq --bed $peak_nopeak_dir/${autosome}.two_consistent_region.nopeak.bed $genome -o ${autosome}.two_consistent_region.nopeak.fa
    done
}

locate_motif() {
    mkdir -p $motif_dir
    cd $motif_dir
    for autosome in $(echo $autosomes | tr " " "\n"); do
        cat $seq_dir/${autosome}.two_consistent_region.peak.fa | $seqkit locate -i -d -f $motif_pattern --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6"\t"$4}' >${autosome}.motif.peak.bed
        cat $seq_dir/${autosome}.two_consistent_region.nopeak.fa | $seqkit locate -i -d -f $motif_pattern --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t"$6"\t"$4}' >${autosome}.motif.nopeak.bed
    done
}
# VARIABLE NAMING TODO:
m6mA_peak_bed="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/peak.bed"
# two_outgroup_consistent_region_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_reichenowi"
two_outgroup_consistent_region_bed_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/two_outgroup_consistent/P_billcollinsi"

peak_nopeak_dir=$two_outgroup_consistent_region_bed_dir/peak_nopeak

seq_dir=$two_outgroup_consistent_region_bed_dir/seq
motif_dir=$two_outgroup_consistent_region_bed_dir/motif_bed
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
motif_pattern="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/across_species/motif"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
get_peak_nopeak_in_two_outgroup_consistent_region
get_seq
locate_motif
