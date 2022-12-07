#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out
#SBATCH --error=%x-%j.err

# FUNCTIONS TODO:
extract_motif_occurrence_zone_seq() {
    cd $motif_occurrence_zone_dir
    pattern=$(basename $motif_occurrence_zone_dir)
    for chrom in $(echo $autosomes | tr " " "\n"); do
        prefix=$(ls | grep $chrom | sed "s/.bed//g")
        $seqkit subseq --bed ${prefix}.bed $genome >${prefix}.fa
    done
    cat *fa >${pattern}.fa
    rm -rf Pf3D7*fa
    cd -
}
# VARIABLE NAMING TODO:
motif_occurrence_zone_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom/relative.0.05"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"

autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"

# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -x
extract_motif_occurrence_zone_seq
