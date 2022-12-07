#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_motif_first_site_appear_in_codon() {
    # sort -k1,1 -k2,2n $cds_each_site_codon_index >${cds_each_site_codon_index}.sort
    motif_each_site=$motif_start_site_dir/nopeak_motif_each_site
    rm -rf $motif_each_site
    for chrom in $(echo $autosomes | tr " " "\n"); do
        awk '$6=="+"{print $1"\t"$2"\t"$2+1};$6=="-"{print $1"\t"$3-1"\t"$3}' $motif_whole_genome_loc/${chrom}.nopeak.genome.motif.bed | sort -k1,1 -k2,2n >$motif_whole_genome_loc/${chrom}.motif.firstSite

        codon_1=$(awk '$3==1' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $motif_whole_genome_loc/${chrom}.motif.firstSite | wc -l)
        codon_2=$(awk '$3==2' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $motif_whole_genome_loc/${chrom}.motif.firstSite | wc -l)
        codon_3=$(awk '$3==3' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $motif_whole_genome_loc/${chrom}.motif.firstSite | wc -l)
        echo $codon_1 $codon_2 $codon_3 >>$motif_each_site
        rm -rf $motif_whole_genome_loc/${chrom}.motif.firstSite
    done
}
# VARIABLE NAMING TODO:
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
motif_whole_genome_loc="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome/motif_GAWGAW_genome_peak_nopeak_loc"
motif_start_site_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/motif_start_site"
cds_each_site_codon_index="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/cds_each_site_codon_index.txt"

bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
get_motif_first_site_appear_in_codon
