#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
working_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/two_outgroup_consistent/single_motif_pattern_GAWGAW_asWhole/GAWGAW/flank_nomatter_strand.remove_close_overlap"
genome="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/genome/PlasmoDB-36_Pfalciparum3D7_Genome.fasta"
chroms="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"

seqkit="/picb/evolgen/users/gushanshan/GenomeAnnotation/seqkit/seqkit"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd $working_dir

for site in $(seq 1 10); do
    cd flank_${site}
    for chrom in $(echo $chroms | tr " " "\n"); do
        $seqkit subseq --bed ${chrom}.peak.bed $genome >${chrom}.peak.fa
    done
    cd ..
done
