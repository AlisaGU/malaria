#!/bin/bash
#SBATCH -n 1
#SBATCH -p hpc
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
get_motif_3_6_site_appear_in_codon123() {
    motif_each_site=$motif_start_site_dir/nopeak_motif_3_6_in_codon123
    rm -rf $motif_each_site
    for chrom in $(echo $autosomes | tr " " "\n"); do

        awk '$6=="+"{print $1"\t"$2+2"\t"$2+3ORS$1"\t"$2+5"\t"$2+6};$6=="-"{print $1"\t"$3-3"\t"$3-2ORS$1"\t"$3-6"\t"$3-5}' $motif_whole_genome_loc/${chrom}.nopeak.genome.motif.bed | sort -k1,1 -k2,2n -u >$motif_whole_genome_loc/${chrom}.nopeak.motif.3_6

        codon1=$(awk '$3==1' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $motif_whole_genome_loc/${chrom}.nopeak.motif.3_6 | wc -l)
        codon1_total=$(awk -v chrom=$chrom '$1==chrom && $3==1' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $chrom_wholeGenome_dir/$chrom.nopeak.bed | wc -l)

        codon2=$(awk '$3==2' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $motif_whole_genome_loc/${chrom}.nopeak.motif.3_6 | wc -l)
        codon2_total=$(awk -v chrom=$chrom '$1==chrom && $3==2' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $chrom_wholeGenome_dir/$chrom.nopeak.bed | wc -l)

        codon3=$(awk '$3==3' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $motif_whole_genome_loc/${chrom}.nopeak.motif.3_6 | wc -l)
        codon3_total=$(awk -v chrom=$chrom '$1==chrom && $3==3' ${cds_each_site_codon_index}.sort | awk '{print $1"\t"$2-1"\t"$2}' | $bedtools intersect -a stdin -b $chrom_wholeGenome_dir/$chrom.nopeak.bed | wc -l)

        echo -e $codon1"\t"$codon1_total"\t"$codon2"\t"$codon2_total"\t"$codon3"\t"$codon3_total >>$motif_each_site

        rm -rf $motif_whole_genome_loc/${chrom}.nopeak.motif.3_6
    done
}
# VARIABLE NAMING TODO:
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
motif_whole_genome_loc="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome/motif_GAWGAW_genome_peak_nopeak_loc"
motif_start_site_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/motif_start_site"
cds_each_site_codon_index="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/ref_genome/pf_3D7/anno/cds_each_site_codon_index.txt"
chrom_wholeGenome_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome"

bedtools=/picb/evolgen/users/gushanshan/software/bedtools/bedtools
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
get_motif_3_6_site_appear_in_codon123
