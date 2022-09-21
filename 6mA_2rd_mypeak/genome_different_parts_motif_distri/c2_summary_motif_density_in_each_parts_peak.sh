#!/bin/bash
#SBATCH -n 1
#SBATCH -A evolgen-grp
#SBATCH -p evolgen
#SBATCH -o %x-%j.out

# FUNCTIONS TODO:
compute_average_motif_count_and_base_count() {
    density_file=$motif_density_dir/peak_motif_density
    rm -rf $density_file

    for chrom in $(echo $autosomes | tr " " "\n"); do
        value=""
        for region in $(echo "intergenic TSS_2KB TTS_2KB genes exons introns" | tr " " "\n"); do
            region_bed=$genome_parts_dir/${region}.bed
            motif_count=$($bedtools intersect -a $motif_whole_genome_loc/${chrom}.peak.genome.motif.bed -b $region_bed | wc -l | awk '{print $1}')
            motif_base_count=$($bedtools intersect -a $motif_whole_genome_loc/${chrom}.peak.genome.motif.bed -b $region_bed | sort -k1,1 -k2,2n | $bedtools merge -i stdin | awk '{print $3-$2}' | awk '{s+=$1}END{print s}')
            # region_len=$(awk -v chrom=$chrom '$1==chrom{print $3-$2}' $region_bed | awk '{s+=$1};END{print s}')
            region_len=$(awk -v chrom=$chrom '$1==chrom' $region_bed | $bedtools intersect -a stdin -b $chrom_wholeGenome_dir/${chrom}.peak.bed | awk '{s+=$3-$2}END{print s}')
            value+=" $(echo $motif_count" "$motif_base_count" "$region_len)"
        done
        echo $value >>$density_file
    done
}
# VARIABLE NAMING TODO:
global_output_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/genome_different_parts_motif_distri"
motif_whole_genome_loc="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome/motif_GAWGAW_genome_peak_nopeak_loc"
genome_parts_dir=$global_output_dir/genome_parts

motif_density_dir=$global_output_dir/motif_density
export autosomes="Pf3D7_01_v3 Pf3D7_02_v3 Pf3D7_03_v3 Pf3D7_04_v3 Pf3D7_05_v3 Pf3D7_06_v3 Pf3D7_07_v3 Pf3D7_08_v3 Pf3D7_09_v3 Pf3D7_10_v3 Pf3D7_11_v3 Pf3D7_12_v3 Pf3D7_13_v3 Pf3D7_14_v3"
export chrom_wholeGenome_dir="/picb/evolgen/users/gushanshan/projects/malaria/dataAndResult/6mA/jiang/2rd/macs2_output/chrom_wholeGenome"

bedtools="/picb/evolgen/users/gushanshan/software/bedtools/bedtools"
# VARIABLE NAMING for test module TODO:

# PROCESS TODO:
# set -ex
mkdir -p $motif_density_dir
cd $motif_density_dir
compute_average_motif_count_and_base_count
