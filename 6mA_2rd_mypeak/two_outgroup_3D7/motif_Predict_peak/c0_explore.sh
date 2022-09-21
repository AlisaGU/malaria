#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:

# VARIABLE NAMING TODO:
autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
set -ex
cd /picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome

file=/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/use_motif_predict_peak/GAWGAW/whole_genome_motif_info.txt
rm -rf $file
for chrom in $(echo $autosomes | tr " " "\n"); do
    peak_len=$(awk '{print $3-$2}' ${chrom}.peak.bed | awk '{s+=$1};END{print s}')
    peak_motif_count=$(wc -l motif_GAWGAW_genome_peak_nopeak_loc/${chrom}.peak.genome.motif.bed | awk '{print $1}')
    peak_merged_motif_base_count=$(sort -k1,1 -k2,2n motif_GAWGAW_genome_peak_nopeak_loc/${chrom}.peak.genome.motif.bed | $bedtools merge -i stdin | awk '{print $3-$2}' | awk '{s+=$1};END{print s}')

    nopeak_len=$(awk '{print $3-$2}' ${chrom}.nopeak.bed | awk '{s+=$1};END{print s}')
    nopeak_motif_count=$(wc -l motif_GAWGAW_genome_peak_nopeak_loc/${chrom}.nopeak.genome.motif.bed | awk '{print $1}')
    nopeak_merged_motif_base_count=$(sort -k1,1 -k2,2n motif_GAWGAW_genome_peak_nopeak_loc/${chrom}.nopeak.genome.motif.bed | $bedtools merge -i stdin | awk '{print $3-$2}' | awk '{s+=$1};END{print s}')
    echo -e $peak_len"\t"$peak_motif_count"\t"$peak_merged_motif_base_count"\t"$nopeak_len"\t"$nopeak_motif_count"\t"$nopeak_merged_motif_base_count >>$file
done
