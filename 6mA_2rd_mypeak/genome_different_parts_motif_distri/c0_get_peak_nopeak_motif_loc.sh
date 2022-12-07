#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_motif() {
    peak_type=$1
    chrom=$2

    $seqkit subseq --quiet --bed $chrom_wholeGenome_dir/${chrom}.${peak_type}.bed $genome | $seqkit locate -i -d -f $motif_pattern --bed | awk '{split($1,coord,/[_-]/);s=coord[4]+$2-1;e=coord[4]+$3-1;print coord[1]"_"coord[2]"_"coord[3]"\t"s"\t"e"\t""GAWGAW""\t""0""\t"$6}' >$motif_GAWGAW_genome_peak_nopeak_loc_dir/${chrom}.${peak_type}.genome.motif.bed
}
export -f get_motif
# VARIABLE NAMING TODO:
export chrom_wholeGenome_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome"
export motif_GAWGAW_genome_peak_nopeak_loc_dir=$chrom_wholeGenome_dir/motif_GAWGAW_genome_peak_nopeak_loc

export genome=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta
export motif_pattern="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/motifs/GAWGAW"

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

export seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
parallel="/home/gushanshan/anaconda3/envs/vscode_r/bin/parallel"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
mkdir -p $motif_GAWGAW_genome_peak_nopeak_loc_dir
$parallel -j 4 get_motif ::: "peak" "nopeak" ::: $autosomes
